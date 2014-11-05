"""
    Created on Apr 11, 2013
    
    This class implements the epitope selection functionality
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
    
    @author: Benjamin Schubert
"""

from __future__ import division

import itertools as itr
import copy

import coopr.environ
from coopr.pyomo import *
from coopr.opt import SolverFactory

from Fred2.Core.Result import EpitopePredictionResult


class OptiTope(object):
    """
        classdocs

            :param EpitopePredictionResult _result: Epitope prediction result object from which the epitope selection should be performed
            :param list(Allele) _alleles: A list of allele object which should be considered during selection and were previously
            used to construct the result object
            :param _k (int): the number of epitopes to select
            :param solver (String): the solver to be used (default glpsol)
    """

    def __init__(self, _results,  threshold=None, k=10, solver="glpsol", verbosity=0):
        """
            Constructor


        """

        #check input data
        if not isinstance(_results, EpitopePredictionResult):
            raise ValueError("first input parameter is not of type EpitopePredictionResult")

        _alleles = copy.deepcopy(_results.columns.values.tolist())

        print map(lambda x: x.locus, _alleles)
        #test if allele prob is set, if not set allele prob uniform
        #if only partly set infer missing values (assuming uniformity of missing value)
        prob = []
        no_prob = []
        for a in _alleles:
            if a.prob is None:
                no_prob.append(a)
            else:
                prob.append(a)

        print no_prob
        if len(no_prob) > 0:
            #group by locus
            no_prob_grouped = {}
            prob_grouped = {}
            for a in no_prob:
                no_prob_grouped.setdefault(a.locus, []).append(a)
            for a in prob:
                prob_grouped.setdefault(a.locus, []).append(a)

            print no_prob_grouped, prob_grouped
            for g, v in no_prob_grouped.iteritems():
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
                print a.name, a.prob

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

        res_df = _results.xs(_results.index.values[0][1], level="Method")
        res_df = res_df[res_df.apply(lambda x: any(x[a] > self.__thresh.get(a.name, -float("inf"))
                                                   for a in res_df.columns), axis=1)]

        for tup in res_df.itertuples():
            p = tup[0]
            seq = str(p)
            peps[seq] = p
            for a, s in itr.izip(res_df.columns, tup[1:]):
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
        for e, v in cons.iteritems():
            try:
                cons[e] = v / total
            except ZeroDivisionError:
                cons[e] = 1
        model = ConcreteModel()

        #set definition
        model.Q = Set(initialize=variations)

        model.E = Set(initialize=set(peps.keys()))

        model.A = Set(initialize=alleles_I.keys())
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
        #model.AntigenCovConst = Constraint(model.Q,
        #                                   rule=lambda model, q: sum(model.x[e] for e in model.E_var[q]) >= model.t_var)
        model.IsAntigenCovConst = Constraint(model.Q,
                                             rule=lambda model, q: sum(model.x[e] for e in model.E_var[q]) >= model.z[q])
        model.MinAntigenCovConst = Constraint(rule=lambda model: sum(model.z[q] for q in model.Q) >= model.t_var)
        model.EpitopeConsConst = Constraint(model.E,
                                            rule=lambda model, e: (1 - model.c[e]) * model.x[e] <= 1 - model.t_c)

        #generate instance
        self.instance = model.create()
        if self.__verbosity > 0:
            print "MODEL INSTANCE"
            self.instance.pprint()

        #deactivate additional constraints and variables
        #params
        self.instance.c.deactivate()
        self.instance.t_c.deactivate()
        self.instance.t_allele.deactivate()
        self.instance.t_var.deactivate()

        #constraints
        #self.instance.IsAlleleCovConst.deactivate()
        self.instance.MinAlleleCovConst.deactivate()
        #self.instance.AntigenCovConst.deactivate()
        #self.instance.IsAntigenCovConst.deactivate()
        self.instance.MinAntigenCovConst.deactivate()
        self.instance.EpitopeConsConst.deactivate()

        #variables
        self.instance.y.deactivate()

    def set_k(self, k):
        """
            sets the number of epitopes to select
            @param k: the number of epitopes
            @type k: int
            @exception OptiTopeException: if the input variable is not in the same domain as the parameter
        """
        tmp = self.instance.k.value
        try:
            getattr(self.instance, str(self.instance.k)).set_value(int(k))
            self.__changed = True
        except:
            self.__changed = False
            getattr(self.instance, str(self.instance.k)).set_value(int(tmp))
            raise Exception('set_k', 'An error has occurred during setting parameter k. Please check if k is integer.')

    def activate_allele_coverage_const(self, minCoverage):
        """
            enables the allele Coverage Constraint

            @param minCoverage (float): percentage of alleles which have to be covered
            @exception EpitopeSelectionException: if the input variable is not in the same domain as the parameter
        """
        # parameter
        mc = self.instance.t_allele.value

        try:
            self.instance.t_allele.activate()
            getattr(self.instance, str(self.instance.t_allele)).set_value(int(len(self.__alleleProb) * minCoverage))
            #variables
            self.instance.y.activate()

            #constraints
            self.instance.IsAlleleCovConst.activate()
            self.instance.MinAlleleCovConst.activate()
            self.__changed = True
        except:
            getattr(self.instance, str(self.instance.t_allele)).set_value(mc)
            self.instance.t_allele.deactivate()
            self.__changed = False
            raise Exception(
                'activate_allele_coverage_const","An error occurred during activation of of the allele coverage constraint. ' +
                'Please check your specified minimum coverage parameter to be in the range of 0.0 and 1.0.')

    def deactivate_allele_coverage_const(self):
        """
            deactivates the allele coverage constraint
        """

        # parameter
        self.__changed = True
        self.instance.t_allele.deactivate()

        #variables
        self.instance.y.deactivate()

        #constraints
        self.instance.IsAlleleCovConst.deactivate()
        self.instance.MinAlleleCovConst.deactivate()

    def activate_antigen_coverage_const(self, t_var):
        """
            activates the variation coverage constraint
            @param t_var: the number of epitopes which have to come from each variation
            @type t_var: int
            @exception EpitopeSelectionException: if the input variable is not in the same domain as the parameter

        """
        tmp = self.instance.t_var.value
        try:
            self.instance.t_var.activate()
            getattr(self.instance, str(self.instance.t_var)).set_value(int(t_var))
            self.instance.IsAntigenCovConst.activate()
            self.instance.MinAntigenCovConst.activate()
            self.__changed = True
        except:
            self.instance.t_var.deactivate()
            getattr(self.instance, str(self.instance.t_var)).set_value(int(tmp))
            self.instance.IsAntigenCovConst.deactivate()
            self.instance.MinAntigenCovConst.deactivate()
            self.__changed = False
            raise Exception("activate_antigen_coverage_const",
                            "An error has occurred during activation of the coverage constraint. Please make sure your input is an integer.")

    def deactivate_antigen_coverage_const(self):
        """
            deactivates the variation coverage constraint
        """
        self.__changed = True
        self.instance.IsAntigenCovConst.deactivate()
        self.instance.MinAntigenCovConst.deactivate()

    def activate_epitope_conservation_const(self, t_c, conservation=None):
        """
            activates the epitope conservation constraint
            @param t_c: the percentage of conservation an epitope has to have.
            @type t_c: float [0.0,1.0]
            :param:dict(Peptide,float) conservation: A dict with key=Peptide specifieying a different conservation score
                                                    for each peptide
            @exception EpitopeSelectionException: if the input variable is not in the same domain as the parameter
        """
        if t_c < 0 or t_c > 1:
            raise Exception("activate_epitope_conservation_const",
                            "The conservation threshold is out of its numerical bound. It has to be between 0.0 and 1.0.")

        self.__changed = True
        self.instance.c.activate()
        self.instance.t_c.activate()
        getattr(self.instance, str(self.instance.t_c)).set_value(float(t_c))
        if conservation is not None:
            for e in self.instance.E:
                if e in conservation:
                    getattr(self.instance, str(self.instance.c))[e] = conservation[e]
                else:
                    getattr(self.instance, str(self.instance.c))[e] = 1.0

        self.instance.EpitopeConsConst.activate()

    def deactivate_epitope_conservation_const(self):
        """
            deactivates epitope conservation constraint
        """
        self.__changed = True
        self.instance.c.deactivate()
        self.instance.t_c.deactivate()
        self.instance.EpitopeConsConst.deactivate()

    def solve(self):
        """
            invokes the selected solver and solves the problem.

            @return returns the optimal epitopes
            @rtype PeptideSet
            @exception EpitopeSelectionException: if the solver raised a problem or the solver is not accessible via the PATH environmental variable.
        """
        if self.__changed:
            try:
                self.instance.x.reset()
                self.instance.y.reset()
                self.instance.preprocess()

                res = self.__solver.solve(self.instance)
                self.instance.load(res)
                if self.__verbosity > 0:
                    res.write(num=1)

                if str(res.Solution.status) != 'optimal':
                    print "Could not solve problem - " + str(res.Solution.status) + ". Please check your settings"
                    sys.exit(-1)

                self.__result = [self.__peptideSet[x] for x in self.instance.x if self.instance.x[x].value == 1.0]
                #self.__result.log_metadata("obj", res.Solution.Objective.Value)

               # DEPRECATED CODE ... Dont know how to give additional information to user
               # if self.__instance.y.active:
               #     self.__result.log_metadata("cov_alleles", AlleleSet(
               #         [Allele(y) for y in self.__instance.y if self.__instance.y[y] == 1.0]))

                self.__changed = False
                return self.__result
            except Exception as e:
                print e
                raise Exception("solve",
                                "An Error has occurred during solving. Please check your settings and if the solver is registered in PATH environment variable.")
        else:
            return self.__result
