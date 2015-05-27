from __future__ import division

import itertools as itr
import warnings
import string
import copy

import coopr.environ
from coopr.pyomo import *
from coopr.opt import SolverFactory
from Fred2.Core.Result import EpitopePredictionResult


class MosaicVaccine(object):
    """
    Implements a mosaic vaccine. The model is based
    on the Orienteering Problem based on this formulation
    http://arxiv.org/pdf/1402.1896.pdf
    couldnt find a better one
    maybe this one here:
    http://users.iems.northwestern.edu/~iravani/Orienteering_IIE.pdf
    or
    http://www.fatih.edu.tr/~jesr/tasgetiren.pdf

    """
    def __init__(self, _results,  threshold=None, k=10, solver="glpsol", verbosity=0):
        """
            Constructor


        """
        def suffixPrefixMatch(x, y, k):
            ''' Return length of longest suffix of x of length at least k that
                matches a prefix of y.  Return 0 if there no suffix/prefix
                match has length at least k.
            '''
            if len(x) < k or len(y) < k:
                return 0
            idx = len(y) # start at the right end of y
            # Search right-to-left in y for length-k suffix of x
            while True:
                hit = string.rfind(y, x[-k:], 0, idx)
                if hit == -1: # not found
                    return len(y)
                ln = hit + k
                # See if match can be extended to include entire prefix of y
                if x[-ln:] == y[:ln]:
                    return len(y)-ln # return length of prefix
                idx = hit + k - 1 # keep searching to left in Y
            return -1


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

        imm = {}
        peps = {}
        alleles_I = {}
        variations = []
        epi_var = {}

        #unstack multiindex df to get normal df based on first prediction method
        #and filter for binding epitopes

        res_df = _results.xs(_results.index.values[0][1], level="Method")
        #print "before",len(res_df)
        res_df = res_df[res_df.apply(lambda x: any(x[a] > self.__thresh.get(a.name, -float("inf"))
                                                   for a in res_df.columns), axis=1)]
        #print "after",len(res_df)


        for tup in res_df.itertuples():
            p = tup[0]
            seq = str(p)
            peps[seq] = p
            for a, s in itr.izip(res_df.columns, tup[1:]):
                if s > self.__thresh.get(a.name, -float("inf")):
                    alleles_I.setdefault(a.name, set()).add(seq)
                imm[seq, a.name] = s

        for a in probs.iterkeys():
            imm['start', a] = 0.0
            imm['stop', a] = 0.0

        self.__peptideSet = peps

        #generate overlapping graph as basis of the problem
        overlapping_graph = {(x,y):suffixPrefixMatch(x,y,1) for x,y in
                             itr.product(peps.iterkeys(),peps.iterkeys()) if x!=y }
        self.overlapping = overlapping_graph
        for k in peps.iterkeys():
            overlapping_graph[("start",k)] = len(k)
            overlapping_graph[(k,"stop")] = 0
        #overlapping_graph[("start","stop")] = 0
        #model here:
        #print overlapping_graph

        ep = peps.keys()
        ep_start = ep + ["start"]
        ep_end = ep +["stop"]
        ep.extend(["start","stop"])

        model = ConcreteModel()

        #Sets
        model.A = Set(initialize=probs.keys())


        model.E = Set(initialize=ep)
        model.E_start = Set(initialize=ep_start)
        model.E_end = Set(initialize=ep_end)
        model.E_prim = Set(initialize=peps.keys())
        model.ExE = Set(initialize=overlapping_graph.iterkeys(),dimen=2)


        #Parameters
        model.p = Param(model.A, initialize=lambda model, a: probs[a])
        model.i = Param(model.E, model.A, initialize=lambda model, e, a: imm[e, a])
        model.LMAX = Param(initialize=self.__k, within=PositiveIntegers, mutable=True)
        model.w_ab = Param(model.ExE, initialize=overlapping_graph)
        model.n = Param(initialize=len(model.E))

        # Variable Definition
        model.x = Var(model.ExE, within=Binary)
        model.u = Var(model.E_end,
                      domain=PositiveIntegers, bounds=(2,model.n))

        # Objective definition
        model.Obj = Objective(
            rule=lambda model: sum(model.x[i,j] * sum( model.i[i, a] for a in model.A)
                                   for i,j in model.ExE),
            sense=maximize)

        #OP specific constraints
        model.StartConnectivity = Constraint(
                                             rule=lambda model: sum(model.x['start',e] for e in model.E_prim) == 1)
        model.EndConnectivity = Constraint(
                                             rule=lambda model: sum(model.x[e,'stop'] for e in model.E_prim) == 1)

        model.Connectivity1 = Constraint(model.E_prim,
                                         rule=lambda model, e: sum(model.x[j,e]
                                                                   for j in model.E_start if j !=e ) <= 1)

        model.Connectivity2 = Constraint(model.E_prim,
                                         rule=lambda model, e: sum(model.x[e,j]
                                                                   for j in model.E_end if j !=e ) <= 1)
        model.TourConstraint = Constraint(model.E_prim,
                                          rule=lambda model, e:sum(model.x[i,e]
                                                                   for i in model.E_start if i !=e ) == sum(model.x[e,j]
                                                                   for j in model.E_end if j !=e ))
        model.LengthConstraint = Constraint(rule=lambda model: sum(model.w_ab[i,j]*model.x[i,j] for i,j in model.ExE)
                                                               <= model.LMAX)

        model.Cardinality = Constraint(filter(lambda x:  x[0] !="start" and x[1] != "start", model.ExE),
                                       rule=lambda model, a, b:
                                                  model.u[a]-model.u[b]+1 <= (model.n -1)*(1-model.x[a, b]))
        #generate instance
        self.instance = model.create()
        if self.__verbosity > 0:
            print "MODEL INSTANCE"
        #    self.instance.pprint()

    def solve(self):
        """
        Solves the Epitope Assembly problem and returns an ordered list of the peptides

        :return: list(Peptide) - An order list of the peptides
        """
        self.instance.x.reset()
        self.instance.u.reset()
        self.instance.preprocess()

        res = self.__solver.solve(self.instance,options="threads=4,mip_cuts_all=2,mpi_cuts_covers=3",)#tee=True)
        self.instance.load(res)
        if self.__verbosity > 0:
            res.write(num=1)

        re = { i:j for i,j in self.instance.ExE if self.instance.x[i,j].value}
        curr = 'start'
        ep = []
        while curr != "stop":
            curr = re[curr]
            if curr == "stop":
                break
            ep.append(curr)
            #print ep
        st = ""
        before = "start"
        for e in ep:
            if before == "start":
                st = e
                before = e
            else:
                st += e[len(e)-self.overlapping[before,e]:]
                before = e
        return ep,st,self.instance.Obj()
        #return [ (i,j) for i,j in self.instance.ExE if self.instance.x[i,j].value]
        #return [ (self.instance.u[u].value, self.__peptideSet[u]) for u in sorted(self.instance.u, key=lambda x: self.instance.u[x].value) if u not in ["start", "stop"]]

    def set_k(self, k):
        """
            sets the number of epitopes to select
            @param k: the number of epitopes
            @type k: int
            @exception OptiTopeException: if the input variable is not in the same domain as the parameter
        """
        tmp = self.instance.LMAX.value
        try:
            getattr(self.instance, str(self.instance.LMAX)).set_value(int(k))
            self.__changed = True
        except:
            self.__changed = False
            getattr(self.instance, str(self.instance.LMAX)).set_value(int(tmp))
            raise Exception('set_k', 'An error has occurred during setting parameter k. Please check if k is integer.')

class MosaicVaccine2(object):
    """
    Implements a mosaic vaccine. The model is based
    on the Orienteering Problem based on this formulation
    http://arxiv.org/pdf/1402.1896.pdf
    couldnt find a better one
    maybe this one here:
    http://users.iems.northwestern.edu/~iravani/Orienteering_IIE.pdf
    or
    http://www.fatih.edu.tr/~jesr/tasgetiren.pdf

    """
    def __init__(self, _results,  threshold=None, k=10, solver="glpsol", verbosity=0):
        """
            Constructor


        """
        def suffixPrefixMatch(x, y, k):
            ''' Return length of longest suffix of x of length at least k that
                matches a prefix of y.  Return 0 if there no suffix/prefix
                match has length at least k.
            '''
            if len(x) < k or len(y) < k:
                return 0
            idx = len(y) # start at the right end of y
            # Search right-to-left in y for length-k suffix of x
            while True:
                hit = string.rfind(y, x[-k:], 0, idx)
                if hit == -1: # not found
                    return len(y)
                ln = hit + k
                # See if match can be extended to include entire prefix of y
                if x[-ln:] == y[:ln]:
                    return len(y)-ln # return length of prefix
                idx = hit + k - 1 # keep searching to left in Y
            return -1


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

        imm = {}
        peps = {}
        alleles_I = {}
        variations = []
        epi_var = {}

        #unstack multiindex df to get normal df based on first prediction method
        #and filter for binding epitopes

        res_df = _results.xs(_results.index.values[0][1], level="Method")
        print "before",len(res_df)
        res_df = res_df[res_df.apply(lambda x: any(x[a] > self.__thresh.get(a.name, -float("inf"))
                                                   for a in res_df.columns), axis=1)]
        print "after",len(res_df)


        for tup in res_df.itertuples():
            p = tup[0]
            seq = str(p)
            peps[seq] = p
            for a, s in itr.izip(res_df.columns, tup[1:]):
                if s > self.__thresh.get(a.name, -float("inf")):
                    alleles_I.setdefault(a.name, set()).add(seq)
                imm[seq, a.name] = s

            prots = set(pr for pr in p.get_all_proteins())
            for prot in prots:
                #variations.append(prot.gene_id)
                epi_var.setdefault(prot.gene_id, set()).add(seq)

        for a in probs.iterkeys():
            imm['start', a] = 0.0
            imm['stop', a] = 0.0

        self.__peptideSet = peps

        #generate overlapping graph as basis of the problem
        overlapping_graph = {(x,y):suffixPrefixMatch(x,y,1) for x,y in
                             itr.product(peps.iterkeys(),peps.iterkeys()) if x!=y }
        self.__overl = overlapping_graph
        for k in peps.iterkeys():
            overlapping_graph[("start",k)] = len(k)
            overlapping_graph[(k,"stop")] = 0
        #overlapping_graph[("start","stop")] = 0
        #model here:
        #print overlapping_graph

        ep = peps.keys()
        ep_start = ep + ["start"]
        ep_end = ep +["stop"]
        ep.extend(["start","stop"])

        model = ConcreteModel()

        #Sets
        model.Q = Set(initialize=set(epi_var.iterkeys()))
        model.A = Set(initialize=probs.keys())

        model.E_var = Set(model.Q, initialize=lambda mode, v: epi_var[v])
        model.A_I = Set(model.A, initialize=lambda model, a: alleles_I[a])

        model.E = Set(initialize=ep)
        model.E_start = Set(initialize=ep_start)
        model.E_end = Set(initialize=ep_end)
        model.E_prim = Set(initialize=peps.keys())
        model.ExE = Set(initialize=overlapping_graph.iterkeys(),dimen=2)


        #Parameters
        model.p = Param(model.A, initialize=lambda model, a: probs[a])
        model.i = Param(model.E, model.A, initialize=lambda model, e, a: imm[e, a])
        model.LMAX = Param(initialize=self.__k, within=PositiveIntegers, mutable=True)
        model.w_ab = Param(model.ExE, initialize=overlapping_graph)
        model.n = Param(initialize=len(model.E))
        model.t_allele = Param(initialize=0, within=NonNegativeIntegers, mutable=True)
        model.t_var = Param(initialize=0, within=NonNegativeIntegers, mutable=True)

        # Variable Definition
        model.x = Var(model.ExE, within=Binary)
        model.y = Var(model.A, within=Binary)
        model.z = Var(model.Q, within=Binary)
        model.u = Var(model.E_end,
                      domain=PositiveIntegers, bounds=(2,model.n))

        # Objective definition
        model.Obj = Objective(
            rule=lambda model: sum(model.x[i,j] * sum( model.i[i, a] for a in model.A)
                                   for i,j in model.ExE),
            sense=maximize)

        #OP specific constraints
        model.StartConnectivity = Constraint(
                                             rule=lambda model: sum(model.x['start',e] for e in model.E_prim) == 1)
        model.EndConnectivity = Constraint(
                                             rule=lambda model: sum(model.x[e,'stop'] for e in model.E_prim) == 1)

        model.Connectivity1 = Constraint(model.E_prim,
                                         rule=lambda model, e: sum(model.x[j,e]
                                                                   for j in model.E_start if j !=e ) <= 1)

        model.Connectivity2 = Constraint(model.E_prim,
                                         rule=lambda model, e: sum(model.x[e,j]
                                                                   for j in model.E_end if j !=e ) <= 1)
        model.TourConstraint = Constraint(model.E_prim,
                                          rule=lambda model, e:sum(model.x[i,e]
                                                                   for i in model.E_start if i !=e ) == sum(model.x[e,j]
                                                                   for j in model.E_end if j !=e ))
        model.LengthConstraint = Constraint(rule=lambda model: sum(model.w_ab[i,j]*model.x[i,j] for i,j in model.ExE)
                                                               <= model.LMAX)

        model.Cardinality = Constraint(filter(lambda x:  x[0] !="start" and x[1] != "start", model.ExE),
                                       rule=lambda model, a, b:
                                                  model.u[a]-model.u[b]+1 <= (model.n -1)*(1-model.x[a, b]))

        #Epitope selection Constraints
        #optional constraints (in basic model they are disabled)
        model.IsAlleleCovConst = Constraint(model.A,
                                            rule=lambda model, a: sum(sum(model.x[e,j] for j in model.E_prim if e != j)
                                                                      for e in model.A_I[a]) >= model.y[a])
        model.MinAlleleCovConst = Constraint(rule=lambda model: sum(model.y[a] for a in model.A) >= model.t_allele)

        model.IsAntigenCovConst = Constraint(model.Q,
                                             rule=lambda model, q: sum(sum(model.x[e,j] for j in model.E_prim if e != j) for e in model.E_var[q]) >= model.z[q])
        model.MinAntigenCovConst = Constraint(rule=lambda model: sum(model.z[q] for q in model.Q) >= model.t_var)

        #generate instance
        self.instance = model.create()
        if self.__verbosity > 0:
            print "MODEL INSTANCE"
        #    self.instance.pprint()

        self.instance.t_allele.deactivate()
        self.instance.t_var.deactivate()

        self.instance.IsAlleleCovConst.deactivate()
        self.instance.MinAlleleCovConst.deactivate()
        self.instance.IsAntigenCovConst.deactivate()
        self.instance.MinAntigenCovConst.deactivate()

        #variables
        self.instance.y.deactivate()
        self.instance.z.deactivate()

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
            self.instance.y.deactivate()

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
            self.instance.z.activate()
            self.__changed = True
        except:
            self.instance.t_var.deactivate()
            getattr(self.instance, str(self.instance.t_var)).set_value(int(tmp))
            self.instance.IsAntigenCovConst.deactivate()
            self.instance.MinAntigenCovConst.deactivate()
            self.instance.z.deactivate()

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
        self.instance.z.deactivate()


    def solve(self):
        """
        Solves the Epitope Assembly problem and returns an ordered list of the peptides

        :return: list(Peptide) - An order list of the peptides
        """
        if self.__changed:
            #try:
                self.instance.x.reset()
                self.instance.y.reset()
                self.instance.preprocess()

                res = self.__solver.solve(self.instance,options="threads=4,mip_cuts_all=2,mip_cuts_covers=3",tee=True)
                self.instance.load(res)
                if self.__verbosity > 0:
                    res.write(num=1)

                if str(res.Solution.status) != 'optimal':
                    print "Could not solve problem - " + str(res.Solution.status) + ". Please check your settings"
                    sys.exit(-1)

                print "covered Antigens ", [q for q in self.instance.Q if self.instance.z[q]]
                re = { i:j for i,j in self.instance.ExE if self.instance.x[i,j].value}
                self.__changed = False
                curr = 'start'
                ep = []
                while curr != "stop":
                    curr = re[curr]
                    if curr == "stop":
                        break
                    ep.append(curr)
                    #print ep
                st = ""
                before = "start"
                for e in ep:
                    if before == "start":
                        st = e
                        before = e
                    else:
                        st += e[len(e)-self.__overl[before,e]:]
                        before = e
                self.__result = (ep,st)
                return ep,st
            #except Exception as e:
            #    print e
            #    raise Exception("solve",
            #                    "An Error has occurred during solving. Please check your settings and if the solver is registered in PATH environment variable.")
        else:
            return self.__result