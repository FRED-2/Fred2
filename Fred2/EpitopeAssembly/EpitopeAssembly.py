# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: EpitopeAssembly.EpitopeAssembly
   :synopsis: This module contains all classes for EpitopeAssembly.
.. moduleauthor:: schubert

"""
from __future__ import division

import os
import subprocess
import warnings
import itertools as itr
import multiprocessing as mp
import copy

from tempfile import NamedTemporaryFile

from pyomo.environ import *
from pyomo.opt import SolverFactory,SolverStatus, TerminationCondition

from Fred2.Core.Base import ACleavageSitePrediction
from Fred2.Core.Peptide import Peptide
from Fred2.CleavagePrediction.PSSM import APSSMCleavageSitePredictor
from Fred2.EpitopePrediction.PSSM import APSSMEpitopePrediction


class EpitopeAssembly(object):
    """
        Implements the epitope assembly approach proposed by Toussaint et al.
        using proteasomal cleavage site prediction and formulating the problem as
        TSP.

        .. note::

            Toussaint, N.C., et al. Universal peptide vaccines - Optimal peptide vaccine design based on viral
            sequence conservation. Vaccine 2011;29(47):8745-8753.

        :param peptides: A list of :class:`~Fred2.Core.Peptide.Peptide` which shell be arranged
        :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`)
        :param pred: A :class:`~Fred2.Core.Base.ACleavageSitePrediction`
        :type pred: :class:`~Fred2.Core.Base.ACleavageSitePredictor`
        :param str solver: Specifies the solver to use (mused by callable by pyomo)
        :param float weight: Specifies how strong unwanted cleavage sites should be punished [0,1], where 0 means they
                             will be ignored, and 1 the sum of all unwanted cleave sites is subtracted from the cleave
                             site between two epitopes
        :param int verbosity: Specifies how verbos the class will be, 0 means normal, >0 debug mode
    """

    def __init__(self, peptides, pred, solver="glpk", weight=0.0, matrix=None, verbosity=0):

        if not isinstance(pred, ACleavageSitePrediction):
            raise ValueError("Cleave site predictor must be of type ACleavageSitePrediction")

        if len(peptides) > 60:
            warnings.warn("The peptide set exceeds 60. Above this level one has to expect " +
                          "considerably long running times due to the complexity of the problem.")

        #Generate model
        #1. Generate peptides for which cleave sites have to be predicted
        #2. generate graph with dummy element
        self.__verbosity = verbosity

        pep_tmp = peptides[:]
        pep_tmp.append("Dummy")
        edge_matrix = {}
        fragments = {}
        seq_to_pep = {}
        self.neo_cleavage = {}
        self.good_cleavage = {}

        if matrix is None:
            for start, stop in itr.combinations(pep_tmp, 2):
                if start == "Dummy" or stop == "Dummy":
                    seq_to_pep[str(start)] = start
                    seq_to_pep[str(stop)] = stop
                    edge_matrix[(str(start), str(stop))] = 0
                    edge_matrix[(str(stop), str(start))] = 0
                else:
                    start_str = str(start)
                    stop_str = str(stop)
                    frag = Peptide(start_str+stop_str)
                    garf = Peptide(stop_str+start_str)

                    fragments[frag] = (start_str, stop_str)
                    fragments[garf] = (stop_str, start_str)

            cleave_pred = pred.predict(fragments.keys())
            #cleave_site_df = cleave_pred.xs((slice(None), (cleavage_pos-1)))
            for i in set(cleave_pred.index.get_level_values(0)):
                fragment = "".join(cleave_pred.ix[i]["Seq"])
                start, stop = fragments[fragment]

                cleav_pos = len(str(start)) - 1
                edge_matrix[(start, stop)] = -1.0 * (cleave_pred.loc[(i, len(str(start)) - 1), pred.name] - weight * sum(
                    cleave_pred.loc[(i, j), pred.name] for j in xrange(cleav_pos-1,cleav_pos+4,1) if j != cleav_pos))

                self.neo_cleavage[(start, stop)] = sum(cleave_pred.loc[(i, j), pred.name] for j in xrange(cleav_pos-1,cleav_pos+4,1) if j != cleav_pos)
                self.good_cleavage[(start, stop)] = cleave_pred.loc[(i, len(str(start)) - 1), pred.name]
        else:
            edge_matrix = matrix
            seq_to_pep = {str(p):p for p in pep_tmp}
            for p in seq_to_pep.iterkeys():
                if p != "Dummy":
                    edge_matrix[(p,"Dummy")] = 0
                    edge_matrix[("Dummy",p)] = 0
        self.__seq_to_pep = seq_to_pep

        #3. initialize ILP
        self.__solver = SolverFactory(solver)
        model = ConcreteModel()

        E = filter(lambda x: x != "Dummy", seq_to_pep.keys())
        model.E = Set(initialize=E)
        model.E_prime = Set(initialize=seq_to_pep.keys())
        model.ExE = Set(initialize=itr.permutations(E,2), dimen=2)

        model.w_ab = Param(model.E_prime, model.E_prime, initialize=edge_matrix)
        model.card = Param(initialize=len(model.E_prime))

        model.x = Var(model.E_prime, model.E_prime, within=Binary)
        model.u = Var(model.E, domain=PositiveIntegers, bounds=(2,model.card))

        model.obj =Objective(
            rule=lambda mode: sum( model.w_ab[a,b]*model.x[a,b] for a in model.E_prime
                                   for b in model.E_prime if a != b),
            sense=minimize)

        model.tour_constraint_1 = Constraint(model.E_prime,
                                             rule=lambda model, a:
                                             sum(model.x[a,b] for b in model.E_prime if a != b) == 1)
        model.tour_constraint_2 = Constraint(model.E_prime,
                                             rule=lambda model, a:
                                             sum(model.x[b,a] for b in model.E_prime if a != b) == 1)
        model.cardinality_constraint = Constraint(model.ExE,
                                                  rule=lambda model, a, b:
                                                  model.u[a]-model.u[b]+1 <= (model.card -1)*(1-model.x[a, b]))

        self.instance = model
        if self.__verbosity > 0:
            print "MODEL INSTANCE"
            self.instance.pprint()

    def solve(self, options=None):
        """
        Solves the Epitope Assembly problem and returns an ordered list of the peptides

        .. note::

            This can take quite long and should not be done for more and 30 epitopes max!

        :param str options: Solver specific options as string (will not be checked for correctness)
        :return: An order list of the :class:`~Fred2.Core.Peptide.Peptide` (based on the string-of-beads ordering)
        :rtype: list(:class:`~Fred2.Core.Peptide.Peptide`)
        """

        options = dict() if options is None else options
        res = self.__solver.solve(self.instance, options=options)
        self.instance.solutions.load_from(res)
        if self.__verbosity > 0:
            res.write(num=1)

        return [self.__seq_to_pep[u] for u in sorted(self.instance.u, key=lambda x: self.instance.u[x].value)]

    def approximate(self):
        """
        Approximates the eptiope assembly problem by applying Lin-Kernighan traveling salesman heuristic

        .. note::

            LKH implementation must be downloaded, compiled, and globally executable.
            Source code can be found here:
            http://www.akira.ruc.dk/~keld/research/LKH/

        :return: An order list of the :class:`~Fred2.Core.Peptide.Peptide` (based on the sting-of-beads ordering)
        :rtype: list(:class:`~Fred2.Core.Peptide.Peptide`)
        """
        tmp_conf = NamedTemporaryFile(delete=False)
        tmp_prob = NamedTemporaryFile(delete=False)
        tmp_out = NamedTemporaryFile(delete=False)


        #write config file:
        tmp_conf.write("PROBLEM_FILE = %s\nOUTPUT_TOUR_FILE = %s\n"%(tmp_prob.name,tmp_out.name))
        tmp_conf.close()
        epis = []
        #write problem file:
        tmp_prob.write("NAME: %s\nTYPE: ATSP\nDIMENSION: %i\nEDGE_WEIGHT_TYPE: EXPLICIT\nEDGE_WEIGHT_FORMAT: FULL_MATRIX\nEDGE_WEIGHT_SECTION\n"%(tmp_prob.name,len(self.instance.E_prime)))
        for i in self.instance.E_prime:
            epis.append(i)
            tmp_prob.write("\t".join("99999999" if i == j else str(int(float(self.instance.w_ab[i,j])*10000)) for j in self.instance.E_prime)+"\n")

        tmp_prob.write("EOF")
        tmp_prob.close()

        #try:
        r = subprocess.call("LKH %s"%tmp_conf.name, shell=True)
        if r == 127:
                raise RuntimeError("LKH is not installed or globally executable.")
        elif r != 0:
                raise RuntimeError("An unknown error occurred for method LKH. "
                                   "Please check whether LKH is globally executable.")
        result = []
        with open(tmp_out.name, "r") as resul_f:
            is_tour = False
            for l in resul_f:
                if is_tour:
                    if int(l.strip()) == -1:
                        break
                    seq = self.__seq_to_pep[epis[int(l.strip())-1]]
                    if seq != "Dummy":
                        result.append(seq)
                elif "TOUR_SECTION" in l:
                    is_tour = True
                else:
                    pass

        tmp_out.close()
        os.remove(tmp_conf.name)
        os.remove(tmp_prob.name)
        os.remove(tmp_out.name)
        return result



########################################################################################################################
def _runs_lexmin(a):
    """
    private used to unpack arguments send to processes
    :param a:
    :return: ei,ej,cleavage_score,imm_score,c1_score,c2_score,non-junction_score
    """
    spacer, cleav, epi, good_cleav, bad_cleav, non_c = _spacer_design(*a)
    return a[0], a[1], cleav, epi, spacer, good_cleav, bad_cleav, non_c


def _spacer_design(ei, ej, k, en, cn, cl_pssm, epi_pssms, cleav_pos, allele_prob, alpha,
                    thresh, solver, beta=0, options=None):
    """
        PRIVATE:
        internal spacer design for a pre-defined spacer length between two epitopes

        :param str ei: start epitope
        :param str ej: end epitope
        :param int k: length of spacer
        :param int en: epitope length
        :param int cn: cleavage-site string length
        :param dict(int:dict(string:float)) cl_pssm: a cleavage site prediction PSSM as dict-of-dicts
        :param dict(int:dict(string:float)) epi_pssm: a epitope prediction PSSM as dict-of-dicts
        :param int cleav_pos: integer specifying at which AA within the epitope of length cn the cleave is predicted
        :param dict(string:float) allele_prob: a dict of HLA alleles as string (i.e. A*02:01) and probabilities [0,1]
        :param float alpha: specifies the first-order influence on the objectives [0,1]
        :param float thresh: specifies at which score a peptide is considered as epitope
        :param string solver: string specifying which ILP solver should be used
        :param dict(str:str) options: solver specific options as keys and parameters as values
        :return: Tuple of ei, ej, spacer (str), cleavage score, immunogenicity score
    """
    options = dict() if options is None else options

    if k <= 0:
        seq = ei+ej
        i = len(ei)-cleav_pos
        g=len(ei)+k-cleav_pos
        c1=sum(cl_pssm[j][seq[i+j]] for j in xrange(cn))+cl_pssm.get(-1,{}).get("con",0)
        c2=sum(cl_pssm[j][seq[g+j]] for j in xrange(cn))+cl_pssm.get(-1,{}).get("con",0)
        non_c = sum(sum(cl_pssm[j][seq[k+j]] for j in xrange(cn) if k != i and k != g)+cl_pssm.get(-1,{}).get("con",0)
                                                        for k in xrange(len(seq)-(cn-1)))

        imm = sum(prob*sum(max(sum(epi_pssms[j,seq[i+j],a] for j in xrange(en))+epi_pssms.get((-1,"con",a),0)-thresh[a],0)
                                                                    for i in xrange(len(seq)-en))
                                                                        for a,prob in allele_prob.iteritems())

        return "",(c1+c2)/2, imm,c1,c2,non_c

    def normalize_pssm(p):
        max_p = -float("inf")
        min_p = float("inf")
        norm = {}
        for i,v in p.iteritems():
            max_tmp = max(v.values())
            min_tmp = min(v.values())
            if max_tmp > max_p:
                max_p = max_tmp
            if min_tmp < min_p:
                min_p = min_tmp
        for i,v in p.iteritems():
            for a,score in v.iteritems():
                norm.setdefault(i,{}).update({a:(score-min_p)/(max_p - min_p)})
        return norm


    cl_pssm_norm = normalize_pssm(cl_pssm)
    alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    model = ConcreteModel()
    le = len(ei)+len(ej)+k
    neg_inf = -float("inf")
    #Sets
    ep = {}
    def sequence_set(model,i):
            if i < len(ei):
                return [ei[i]]
            elif i < len(ei)+k:
                return alphabet
            else:
                return [ej[i-len(ei)-k]]
    model.A = Set(initialize=allele_prob.keys())
    model.C = Set(initialize=range(cn))
    model.EN = Set(initialize=range(en))
    model.L = Set(initialize=range(le))
    model.Sigma = Set(initialize=alphabet)
    model.S = Set(model.L, initialize=sequence_set)
    model.AUX = Set(dimen=2, initialize=lambda model: [(i,a) for i in xrange(le) for a in model.S[i]])
    model.R = Set(initialize=range(le-(en-1)))

    #param
    model.f = Param(model.C, model.Sigma, initialize=lambda model,i,a: cl_pssm_norm[i][a])
    model.ci = Param(initialize=len(ei)-cleav_pos)
    model.cj = Param(initialize=len(ei)+k-cleav_pos)
    model.p = Param(model.A, initialize=lambda model, m: allele_prob[m])
    model.bi = Param(model.A, initialize=lambda model, m: epi_pssms.get((-1,"con",m),0))
    model.bc = Param(initialize=cl_pssm_norm.get(-1,{}).get("con",0))
    #epitope part
    model.i = Param(model.EN, model.Sigma,model.A, initialize=lambda model,i,a,m: epi_pssms[i,a,m])
    model.tau_epi = Param(initialize=10**6,mutable=True)
    model.tau_cleav = Param(initialize=-10**6, mutable=True)
    model.t_a = Param(model.A, initialize=lambda model, a: thresh.get(a, 0))

    #Variables
    model.x = Var(model.AUX,domain=Binary)
    model.y = Var(model.R,model.A,domain=NonNegativeReals)

    #objective linear
    model.obj_cleav = Objective(rule=lambda model: 0.5*(sum( model.f[i,a]*model.x[model.ci+i,a] for i in model.C for a in model.S[model.ci+i] )
                              + sum(model.f[j,a]*model.x[model.cj+j,a] for j in model.C for a in model.S[model.cj+j])+2*model.bc),sense=maximize)

    model.obj_epi = Objective(rule=lambda model: sum(model.y[i,a]*model.p[a] for a in model.A
                                                     for i in model.R), sense=minimize)

    model.obj_non_cleav = Objective(rule=lambda model: sum( model.f[j,a]*model.x[j+i,a] for i in xrange(le-(cn-1))
                                                            for j in model.C
                                                                for a in model.S[i+j]
                                                                    if i != model.ci and i != model.cj)+(le-(cn-1)-2)*model.bc, sense=minimize)



    #constraints
    model.cons = Constraint(model.L,rule=lambda model, i: sum(model.x[i,a] for a in model.S[i]) == 1)

    model.max_imm_c = Constraint(model.R,model.A,rule=lambda model, i, m:
                                            model.y[i,m] >= sum( model.x[i+j,a]*model.i[j,a,m]
                                                                for j in model.EN for a in model.S[i+j])+model.bi[m]-model.t_a[m])

    ##neo-epitope constraint
    model.c_epi = Constraint(rule=lambda model:sum(model.y[i,a]*model.p[a] for a in model.A
                                                   for i in model.R) <= model.tau_epi)

    #cleavage constraint
    model.c_cleavage = Constraint(rule=lambda model: 0.5*(sum( model.f[i,a]*model.x[model.ci+i,a] for i in model.C for a in model.S[model.ci+i] )
                              + sum(model.f[j,a]*model.x[model.cj+j,a] for j in model.C for a in model.S[model.cj+j])+2*model.bc) >= model.tau_cleav)

    instance = model
    solver = SolverFactory(solver)

    instance.obj_epi.deactivate()
    instance.obj_non_cleav.deactivate()
    instance.c_epi.deactivate()
    instance.c_cleavage.deactivate()

    res = solver.solve(instance, options=options)#, tee=True)

    if (res.solver.status == SolverStatus.ok) and (res.solver.termination_condition == TerminationCondition.optimal):
        instance.solutions.load_from(res)
        #print "In Second objective ",options
        obj_cleav = instance.obj_cleav()

        instance.obj_cleav.deactivate()
        instance.obj_epi.activate()
        instance.c_cleavage.activate()

        #set bound of now inactive objective
        getattr(instance, "tau_cleav").set_value(alpha*obj_cleav)
        #instance.pprint()

        res2 = solver.solve(instance, options=options)#, tee=True)
        if (res2.solver.status == SolverStatus.ok) and (
            res2.solver.termination_condition == TerminationCondition.optimal):
            instance.solutions.load_from(res2)

            if beta:

                #print "In thrid objective"
                obj_imm = instance.obj_epi()

                instance.obj_epi.deactivate()
                instance.obj_non_cleav.activate()
                instance.c_epi.activate()

                getattr(instance, "tau_epi").set_value((2-beta)*obj_imm)
                #print "imm",obj_imm,"tau_epi", (2-beta)*obj_imm

                res3 = solver.solve(instance, options=options)#, tee=True)
                if (res3.solver.status == SolverStatus.ok) and (res3.solver.termination_condition == TerminationCondition.optimal):
                    instance.load(res3)
                    ci = float(sum(cl_pssm[i][a]*instance.x[model.ci+i,a] for i in instance.C for a in instance.S[instance.ci+i] ))+cl_pssm.get(-1,{}).get("con",0)
                    cj = float(sum(cl_pssm[j][a]*instance.x[model.cj+j,a] for j in instance.C for a in instance.S[instance.cj+j]))+cl_pssm.get(-1,{}).get("con",0)
                    imm = float(sum(instance.y[i,a]*instance.p[a] for a in instance.A for i in instance.R))
                    non_c = float(sum(cl_pssm[j][a]*instance.x[j+i,a] for i in xrange(le-(cn-1))
                                                            for j in instance.C
                                                                for a in instance.S[i+j]
                                                                    if i != instance.ci and i != instance.cj))

                    return "".join([a for i in xrange(len(ei), len(ei) + k) for a in instance.S[i] if
                            instance.x[i, a].value]), float(ci+cj)/2, imm, float(ci),float(cj),non_c
                else:
                    raise RuntimeError("Problem could not be solved. Please check your input.")
            else:
                ci = float(sum(cl_pssm[i][a]*instance.x[model.ci+i,a].value for i in instance.C for a in instance.S[instance.ci+i] ))+cl_pssm.get(-1,{}).get("con",0)
                cj = float(sum(cl_pssm[j][a]*instance.x[model.cj+j,a].value for j in instance.C for a in instance.S[instance.cj+j]))+cl_pssm.get(-1,{}).get("con",0)
                non_c = float(sum(cl_pssm[j][a]*instance.x[j+i,a].value for i in xrange(le-(cn-1))
                                                            for j in instance.C
                                                                for a in instance.S[i+j]
                                                                    if i != instance.ci and i != instance.cj))
                return "".join([a for i in xrange(len(ei), len(ei) + k) for a in instance.S[i] if
                            instance.x[i, a].value]), float(ci+cj)/2, instance.obj_epi(),float(ci),float(cj),non_c
        else:
            raise RuntimeError("Problem could not be solved. Please check your input.")
    else:
        raise RuntimeError("Problem could not be solved. Please check your input.")


class EpitopeAssemblyWithSpacer(object):
    """

        Implements the epitope assembly approach proposed by Toussaint et al.
        using proteasomal cleavage site prediction and formulating the problem as
        TSP.

        It also extends it by optimal spacer design.
        (currently only allowed with PSSM cleavage site and epitope prediction)

        The ILP model is implemented. So be reasonable with the size of epitope to be arranged.
    """

    def __init__(self, peptides, cleav_pred, epi_pred, alleles, k=5, en=9, threshold=None, solver="glpk", alpha=0.99,
                 beta=0, verbosity=0):
        """

        :param peptides: A list of :class:`~Fred2.Core.Peptide.Peptide` which shell be arranged
        :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`)
        :param cleav_pred: A :class:`~Fred2.CleavagePrediction.PSSM.APSSMCleavageSitePredictor` (PSSM only)
        :type cleav_pred: :class:`~Fred2.Core.Base.ACleavageSitePredictor`
        :param epi_pred: A :class:`~Fred2.EpitopePrediction.PSSM.APSSMEpitopePrediction` (PSSM only)
        :type epi_pred: :class:`~Fred2.Core.Base.AEpitopePredictor`
        :param alleles: A list of :class:`~Fred2.Core.Allele.Allele` for which predictions should be made
        :type alleles: list(:class:`~Fred2.Core.Allele.Allele`)
        :param int k: The maximal length of a spacer
        :param int en: Length of epitopes
        :param dict(str,float) threshold: A dictionary specifying the epitope prediction threshold for each
                                          :class:`~Fred2.Core.Allele.Allele`
        :param str solver: Specifies the solver to use (must be callable by pyomo)
        :param float alpha: Specifies how how much junction-cleavage score can be sacrificed  to gain lower
                            neo-immunogenicity
        :param float beta: Specifies how how much noe-immunogenicity score can be sacrificed to gain lower non-junction
                           cleavage score
        :param int verbosity: Specifies how verbos the class will be, 0 means normal, >0 debug mode
        """

        #test input
        if not isinstance(cleav_pred, APSSMCleavageSitePredictor):
            raise ValueError("Second input must be a PSSM-based cleavage site predictor.")

        if not isinstance(epi_pred, APSSMEpitopePrediction):
            raise ValueError("Third input must be a PSSM-based epitope predictor.")

        if en not in epi_pred.supportedLength:
            raise ValueError("Specified epitope length of en=%i is not supported by %s"%(en,epi_pred.name))

        _alleles = [copy.deepcopy(a) for a in alleles if a.name in epi_pred.supportedAlleles]

        if not _alleles:
            raise ValueError("Specified alleles are not supported by %s"%epi_pred.name)

        #infere probability if not already set

        #test if allele prob is set, if not set allele prob uniform
        #if only partly set infer missing values (assuming uniformity of missing value)
        prob = []
        no_prob = []
        for a in _alleles:
            if a.prob is None:
                no_prob.append(a)
            else:
                prob.append(a)

        #print no_prob
        if len(no_prob) > 0:
            #group by locus
            no_prob_grouped = {}
            prob_grouped = {}
            for a in no_prob:
                no_prob_grouped.setdefault(a.locus, []).append(a)
            for a in prob:
                prob_grouped.setdefault(a.locus, []).append(a)

            #print no_prob_grouped, prob_grouped
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


        self.spacer = {}

        #start constructing model
        self.__solver = solver
        self.__verbosity = verbosity
        self.__changed = True
        self.__k = k
        self.__result = None
        self.__thresh = {a.name:0 for a in alleles} if threshold is None else threshold
        self.__alleles = _alleles
        self.__epi_pred = epi_pred
        self.__clev_pred = cleav_pred
        self.__en = en
        self.__alpha = alpha
        self.__beta = beta
        self.__peptides = peptides
        #model construction for spacer design

    def solve(self, start=0, threads=None, options=None):
        """
        Solve the epitope assembly problem with spacers optimally using integer linear programming.

        .. note::

            This can take quite long and should not be done for more and 30 epitopes max!
            Also, one has to disable pre-solving steps in order to use this model.

        :param int start: Start length for spacers (default 0).
        :param int threads: Number of threads used for spacer design.
                            Be careful, if options contain solver threads it will allocate threads*solver_threads cores!
        :param dict(str, str) options: Solver specific options as keys and parameters as values
        :return: A list of ordered :class:`~Fred2.Core.Peptide.Peptide`
        :rtype: list(:class:`~Fred2.Core.Peptide.Peptide`)
        """
        def __load_model(name, model):
            return getattr(__import__("Fred2.Data.pssms."+name+".mat."+model, fromlist=[model]), model)

        options = dict() if options is None else options
        threads = mp.cpu_count() if threads is None else threads
        pool = mp.Pool(threads)


        #prepare parameters
        cn = min(self.__clev_pred.supportedLength)
        cl_pssm = __load_model(self.__clev_pred.name, self.__clev_pred.name+"_"+str(cn))
        cleav_pos = self.__clev_pred.cleavagePos
        en = self.__en
        epi_pssms = {}
        allele_prob = {}
        for a in self.__alleles:
            allele_prob[a.name] = a.prob
            pssm = __load_model(self.__epi_pred.name, "%s_%i"%(self.__epi_pred.convert_alleles([a])[0], en))
            for j, v in pssm.iteritems():
                for aa, score in v.iteritems():
                    if self.__epi_pred.name in ["SMM", "SMMPMBEC"]:
                        epi_pssms[j, aa, a.name] = 1/10. - math.log(math.pow(10, score), 50000)
                        self.__thresh = {k: 1-math.log(v, 50000) for k, v in self.__thresh.iteritems()}
                    else:
                        epi_pssms[j, aa, a.name] = score

        #print "run spacer designs in parallel using multiprocessing"
        res = map(_runs_lexmin, ((str(ei), str(ej), i, en, cn, cl_pssm, epi_pssms, cleav_pos, allele_prob,
                                      self.__alpha, self.__thresh, self.__solver, self.__beta, options)
                                      for i in xrange(start, self.__k+1)
                                      for ei, ej in itr.product(self.__peptides, repeat=2) if ei != ej))
        pool.close()
        pool.join()


        opt_spacer = {}
        adj_matrix = {}
        inf = float("inf")
        #print res
        #print "find best scoring spacer for each epitope pair"
        for ei, ej, score, epi, spacer, c1, c2, non_c in res:
                #print ei,spacer,ej,min(c1,c2),c1,c2
                if adj_matrix.get((ei, ej), inf) > -min(c1,c2):
                    adj_matrix[(ei, ej)] = -min(c1,c2)
                    opt_spacer[(ei, ej)] = spacer

        self.spacer = opt_spacer
        #print "solve assembly with generated adjacency matrix"
        assembler = EpitopeAssembly(self.__peptides, self.__clev_pred, solver=self.__solver, matrix=adj_matrix)
        res = assembler.solve(options=options)

        #generate output
        sob = []
        for i in xrange(len(res)-1):
            ei = str(res[i])
            ej = str(res[i+1])
            if not i:
                sob.append(Peptide(ei))
            sob.append(Peptide(opt_spacer[ei,ej]))
            sob.append(Peptide(ej))
        return sob

    def approximate(self, start=0, threads=1, options=None):
        """
        Approximates the Eptiope Assembly problem by applying Lin-Kernighan traveling salesman heuristic

        LKH implementation must be downloaded, compiled, and globally executable.

        Source code can be found here:
        http://www.akira.ruc.dk/~keld/research/LKH/

        :param int start: Start length for spacers (default 0).
        :param int threads: Number of threads used for spacer design. Be careful, if options contain solver threads it
                            will allocate threads*solver_threads cores!
        :param dict(str, str) options: Solver specific options (threads for example)
        :return: A list of ordered :class:`~Fred2.Core.Peptide.Peptide`
        :rtype: list(:class:`~Fred2.Core.Peptide.Peptide`)
        """
        def __load_model(data, name, length):
            model = "%s_%s"%(name, str(length))
            return getattr(__import__(data, fromlist=[model]), model)

        options = dict() if options is None else options
        threads = mp.cpu_count() if threads is None else threads
        pool = mp.Pool(threads)


        #prepare parameters
        cn = min(self.__clev_pred.supportedLength)
        cl_pssm = __load_model("Fred2.Data.CleaveagePSSMMatrices", self.__clev_pred.name,cn)
        cleav_pos = self.__clev_pred.cleavagePos
        en = self.__en
        epi_pssms = {}
        allele_prob = {}
        for a in self.__alleles:

            try:
                pssm = __load_model("Fred2.Data.EpitopePSSMMatrices",
                                    self.__epi_pred.name, "%s_%i"%(self.__epi_pred.convert_alleles([a])[0], en))
                allele_prob[a.name] = a.prob
                for j, v in pssm.iteritems():
                    for aa, score in v.iteritems():
                        epi_pssms[j, aa, a.name] = score
            except AttributeError:
                #del self.__thresh[a.name]
                continue

        if not epi_pssms:
            raise ValueError("Selected alleles with epitope length are not supported by the prediction method.")

        #print "run spacer designs in parallel using multiprocessing"
        res = pool.map(_runs_lexmin, ((str(ei), str(ej), i, en, cn, cl_pssm, epi_pssms, cleav_pos, allele_prob,
                                       self.__alpha, self.__thresh, self.__solver, self.__beta, options)
                                      for i in xrange(start, self.__k+1)
                                      for ei, ej in itr.product(self.__peptides, repeat=2) if ei != ej))
        pool.close()
        pool.join()


        opt_spacer = {}
        adj_matrix = {}
        inf = float("inf")
        #print res
        #print "find best scoring spacer for each epitope pair"
        for ei, ej, score, epi, spacer, c1, c2, non_c in res:
                if adj_matrix.get((ei, ej), inf) > -min(c1, c2):
                    adj_matrix[(ei, ej)] = -min(c1, c2)
                    opt_spacer[(ei, ej)] = spacer

        self.spacer = opt_spacer
        #print "solve assembly with generated adjacency matrix"
        assembler = EpitopeAssembly(self.__peptides, self.__clev_pred, solver=self.__solver, matrix=adj_matrix)
        res = assembler.approximate()

        #generate output
        sob = []
        for i in xrange(len(res)-1):
            ei = str(res[i])
            ej = str(res[i+1])
            if not i:
                sob.append(Peptide(ei))
            sob.append(Peptide(opt_spacer[ei,ej]))
            sob.append(Peptide(ej))
        return sob