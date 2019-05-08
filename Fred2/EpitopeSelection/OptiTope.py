# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
   .. module:: OptiTope
   :synopsis:  This class implements the epitope selection functionality
    of OptiTope published by Toussaint et al. [1].
    
    This module builds upon Coopr's Pyomo, an embedded algebraic modeling
    languages [2].
    
    It allows to (de)select specific constraints of the underlying
    ILP and to solve the specific problem with a MIP solver of choice
    
    
    [1] N. C. Toussaint and O. Kohlbacher. OptiTope--a web server for the selection of
    an optimal set of peptides for epitope-based vaccines. Nucleic Acids Res,
    2009, 37, W617-W622
    [2] Pyomo - Optimization Modeling in Python. William E. Hart, Carl Laird,
    Jean-Paul Watson and David L. Woodruff. Springer, 2012.
    
.. moduleauthor:: schubert

"""




import itertools as itr
import copy
import math

from pyomo.environ import ConcreteModel, Set, Param, Var, Constraint, PositiveIntegers, \
                          Binary, NonNegativeIntegers, Objective, maximize, NonNegativeReals
from pyomo.opt import SolverFactory, TerminationCondition

from Fred2.Core.Result import EpitopePredictionResult


class OptiTope(object):
    """
    This class implements the epitope selection functionality
    of OptiTope published by Toussaint et al. [1].

    This module builds upon Pyomo, an embedded algebraic modeling
    languages [2].

    It allows to (de)select specific constraints of the underlying
    ILP and to solve the specific problem with a MIP solver of choice

    .. note::

        [1] N. C. Toussaint and O. Kohlbacher. OptiTope--a web server for the selection of
        an optimal set of peptides for epitope-based vaccines. Nucleic Acids Res,
        2009, 37, W617-W622
        [2] Pyomo - Optimization Modeling in Python. William E. Hart, Carl Laird,
        Jean-Paul Watson and David L. Woodruff. Springer, 2012.
    """

    def __init__(self, results,  threshold=None, k=10, solver="glpk", verbosity=0):
        """
        :param result: Epitope prediction result object from which the epitope selection should be performed
        :type result: :class:`~Fred2.Core.Result.EpitopePredictionResult`
        :param dict(str,float) threshold: A dictionary scoring the binding thresholds for each HLA
                                          :class:`~Fred2.Core.Allele.Allele` key = allele name; value = the threshold
        :param int k: The number of epitopes to select
        :param str solver: The solver to be used (default glpk)
        :param int verbosity: Integer defining whether additional debugg prints are made >0 => debug mode
        """

        #check input data
        if not isinstance(results, EpitopePredictionResult):
            raise ValueError("first input parameter is not of type EpitopePredictionResult")

        _alleles = copy.deepcopy(results.columns.values.tolist())

        #test if allele prob is set, if not set allele prob uniform
        #if only partly set infer missing values (assuming uniformity of missing value)
        prob = []
        no_prob = []
        for a in _alleles:
            if a.prob is None:
                no_prob.append(a)
            else:
                prob.append(a)

        if len(no_prob) > 0:
            #group by locus
            no_prob_grouped = {}
            prob_grouped = {}
            for a in no_prob:
                no_prob_grouped.setdefault(a.locus, []).append(a)
            for a in prob:
                prob_grouped.setdefault(a.locus, []).append(a)

            for g, v in no_prob_grouped.items():
                total_loc_a = len(v)
                if g in prob_grouped:
                    remaining_mass = 1.0 - sum(a.prob for a in prob_grouped[g])
                    for a in v:
                        a.prob = remaining_mass/total_loc_a
                else:
                    for a in v:
                        a.prob = 1.0/total_loc_a
        probs = {a.name:a.prob for a in _alleles}
        if verbosity:
            for a in _alleles:
                print(a.name, a.prob)

        #start constructing model
        self.__solver = SolverFactory(solver)
        self.__verbosity = verbosity
        self.__changed = True
        self.__alleleProb = _alleles
        self.__k = k
        self.__result = None
        self.__thresh = {} if threshold is None else threshold

        # Variable, Set and Parameter preparation
        alleles_I = {}
        variations = []
        epi_var = {}
        imm = {}
        peps = {}
        cons = {}

        #unstack multiindex df to get normal df based on first prediction method
        #and filter for binding epitopes
        method = results.index.values[0][1]
        res_df = results.xs(results.index.values[0][1], level="Method")
        res_df = res_df[res_df.apply(lambda x: any(x[a] > self.__thresh.get(a.name, -float("inf"))
                                                   for a in res_df.columns), axis=1)]

        for tup in res_df.itertuples():
            p = tup[0]
            seq = str(p)
            peps[seq] = p
            for a, s in itr.izip(res_df.columns, tup[1:]):
                if method in ["smm", "smmpmbec", "arb", "comblibsidney"]:
                    try:
                        thr = min(1., max(0.0, 1.0 - math.log(self.__thresh.get(a.name),
                                                      50000))) if a.name in self.__thresh else -float("inf")
                    except:
                        thr = 0

                    if s >= thr:
                        alleles_I.setdefault(a.name, set()).add(seq)
                    imm[seq, a.name] = min(1., max(0.0, 1.0 - math.log(s, 50000)))
                else:
                    if s > self.__thresh.get(a.name, -float("inf")):
                        alleles_I.setdefault(a.name, set()).add(seq)
                    imm[seq, a.name] = s

            prots = set(pr for pr in p.get_all_proteins())
            cons[seq] = len(prots)
            for prot in prots:
                variations.append(prot.gene_id)
                epi_var.setdefault(prot.gene_id, set()).add(seq)
        self.__peptideSet = peps

        #calculate conservation
        variations = set(variations)
        total = len(variations)
        for e, v in cons.items():
            try:
                cons[e] = v / total
            except ZeroDivisionError:
                cons[e] = 1
        model = ConcreteModel()

        #set definition
        model.Q = Set(initialize=variations)

        model.E = Set(initialize=set(peps.keys()))

        model.A = Set(initialize=list(alleles_I.keys()))
        model.E_var = Set(model.Q, initialize=lambda mode, v: epi_var[v])
        model.A_I = Set(model.A, initialize=lambda model, a: alleles_I[a])


        #parameter definition
        model.k = Param(initialize=self.__k, within=PositiveIntegers, mutable=True)
        model.p = Param(model.A, initialize=lambda model, a: probs[a])

        model.c = Param(model.E, initialize=lambda model, e: cons[e],mutable=True)

        #threshold parameters
        model.i = Param(model.E, model.A, initialize=lambda model, e, a: imm[e, a])
        model.t_allele = Param(initialize=0, within=NonNegativeIntegers, mutable=True)
        model.t_var = Param(initialize=0, within=NonNegativeIntegers, mutable=True)
        model.t_c = Param(initialize=0.0, within=NonNegativeReals, mutable=True)

        # Variable Definition
        model.x = Var(model.E, within=Binary)
        model.y = Var(model.A, within=Binary)
        model.z = Var(model.Q, within=Binary)

        # Objective definition
        model.Obj = Objective(
            rule=lambda model: sum(model.x[e] * sum(model.p[a] * model.i[e, a] for a in model.A) for e in model.E),
            sense=maximize)


        #Obligatory Constraint (number of selected epitopes)
        model.NofSelectedEpitopesCov = Constraint(rule=lambda model: sum(model.x[e] for e in model.E) <= model.k)

        #optional constraints (in basic model they are disabled)
        model.IsAlleleCovConst = Constraint(model.A,
                                            rule=lambda model, a: sum(model.x[e] for e in model.A_I[a]) >= model.y[a])
        model.MinAlleleCovConst = Constraint(rule=lambda model: sum(model.y[a] for a in model.A) >= model.t_allele)

        model.IsAntigenCovConst = Constraint(model.Q,
                                             rule=lambda model, q: sum(model.x[e] for e in model.E_var[q]) >= model.z[q])
        model.MinAntigenCovConst = Constraint(rule=lambda model: sum(model.z[q] for q in model.Q) >= model.t_var)
        model.EpitopeConsConst = Constraint(model.E,
                                            rule=lambda model, e: (1 - model.c[e]) * model.x[e] <= 1 - model.t_c)

        #generate instance
        self.instance = model
        if self.__verbosity > 0:
            print("MODEL INSTANCE")
            self.instance.pprint()

        #constraints
        self.instance.IsAlleleCovConst.deactivate()
        self.instance.MinAlleleCovConst.deactivate()
        self.instance.IsAntigenCovConst.deactivate()
        self.instance.MinAntigenCovConst.deactivate()
        self.instance.EpitopeConsConst.deactivate()


    def set_k(self, k):
        """
            Sets the number of epitopes to select

            :param int k: The number of epitopes
            :raises ValueError: If the input variable is not in the same domain as the parameter
        """
        tmp = self.instance.k.value
        try:
            getattr(self.instance, str(self.instance.k)).set_value(int(k))
            self.__changed = True
        except ValueError:
            self.__changed = False
            getattr(self.instance, str(self.instance.k)).set_value(int(tmp))
            raise ValueError('set_k', 'An error has occurred during setting parameter k. Please check if k is integer.')

    def activate_allele_coverage_const(self, minCoverage):
        """
            Enables the allele coverage constraint

            :param float minCoverage: Percentage of alleles which have to be covered [0,1]
            :raises ValueError: If the input variable is not in the same domain as the parameter
        """
        # parameter
        mc = self.instance.t_allele.value

        try:
            getattr(self.instance, str(self.instance.t_allele)).set_value(int(len(self.__alleleProb) * minCoverage))
            #variables

            #constraints
            self.instance.IsAlleleCovConst.activate()
            self.instance.MinAlleleCovConst.activate()
            self.__changed = True
        except ValueError:
            getattr(self.instance, str(self.instance.t_allele)).set_value(mc)
            self.__changed = False
            raise ValueError(
                'activate_allele_coverage_const","An error occurred during activation of of the allele coverage constraint. ' +
                'Please check your specified minimum coverage parameter to be in the range of 0.0 and 1.0.')

    def deactivate_allele_coverage_const(self):
        """
            Deactivates the allele coverage constraint
        """

        # parameter
        self.__changed = True

        #constraints
        self.instance.IsAlleleCovConst.deactivate()
        self.instance.MinAlleleCovConst.deactivate()

    def activate_antigen_coverage_const(self, t_var):
        """
            Activates the variation coverage constraint

            :param int t_var: The number of epitopes which have to come from each variation
            :raises ValueError: If the input variable is not in the same domain as the parameter

        """
        tmp = self.instance.t_var.value
        try:
            self.instance.z
            getattr(self.instance, str(self.instance.t_var)).set_value(int(len(self.instance.Q)*t_var))
            self.instance.IsAntigenCovConst.activate()
            self.instance.MinAntigenCovConst.activate()
            self.__changed = True
        except ValueError:
            getattr(self.instance, str(self.instance.t_var)).set_value(int(tmp))
            self.instance.IsAntigenCovConst.deactivate()
            self.instance.MinAntigenCovConst.deactivate()
            self.__changed = False
            raise ValueError("activate_antigen_coverage_const",
                            "An error has occurred during activation of the coverage constraint. Please make sure your input is an integer.")

    def deactivate_antigen_coverage_const(self):
        """
            Deactivates the variation coverage constraint
        """
        self.__changed = True
        self.instance.IsAntigenCovConst.deactivate()
        self.instance.MinAntigenCovConst.deactivate()

    def activate_epitope_conservation_const(self, t_c, conservation=None):
        """
            Activates the epitope conservation constraint

            :param float t_c: The percentage of conservation an epitope has to have [0.0,1.0].
            :param: conservation: A dict with key=:class:`~Fred2.Core.Peptide.Peptide` specifying a different
                                  conservation score for each :class:`~Fred2.Core.Peptide.Peptide`
            :type conservation: dict(:class:`~Fred2.Core.Peptide.Peptide`,float)
            :raises ValueError: If the input variable is not in the same domain as the parameter
        """
        if t_c < 0 or t_c > 1:
            raise ValueError("activate_epitope_conservation_const",
                            "The conservation threshold is out of its numerical bound. It has to be between 0.0 and 1.0.")

        self.__changed = True
        getattr(self.instance, str(self.instance.t_c)).set_value(float(t_c))
        if conservation is not None:
            for e in self.instance.E:
                if e in conservation:
                    getattr(self.instance, str(self.instance.c))[e] = conservation[e]
                else:
                    getattr(self.instance, str(self.instance.c))[e] = 0.0

        self.instance.EpitopeConsConst.activate()

    def deactivate_epitope_conservation_const(self):
        """
            Deactivates epitope conservation constraint
        """
        self.__changed = True
        self.instance.EpitopeConsConst.deactivate()

    def solve(self, options=None):
        """
            Invokes the selected solver and solves the problem

            :param dict(str,str) options: A dictionary of solver specific options as keys and their parameters as values
            :return Returns the optimal epitopes as list of :class:`~Fred2.Core.Peptide.Peptide` objectives
            :rtype: list(:class:`~Fred2.Core.Peptide.Peptide`)
            :raise RuntimeError: If the solver raised a problem or the solver is not accessible via the PATH
                                 environmental variable.
        """
        options = dict() if options is None else options

        if self.__changed:
            try:

                res = self.__solver.solve(self.instance, options=options)
                self.instance.solutions.load_from(res)
                if self.__verbosity > 0:
                    res.write(num=1)

                if res.solver.termination_condition != TerminationCondition.optimal:
                    raise RuntimeError("Could not solve problem - " + str(res.Solution.status) + ". Please check your settings")

                self.__result = [self.__peptideSet[x] for x in self.instance.x if self.instance.x[x].value == 1.0]
                #self.__result.log_metadata("obj", res.Solution.Objective.Value)

                self.__changed = False
                return self.__result
            except Exception as e:
                raise RuntimeError("solve",
                                "An Error has occurred during solving. Please check your settings and if the solver is registered in PATH environment variable.")
        else:
            return self.__result
