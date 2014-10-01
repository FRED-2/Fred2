# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
from __future__ import division
__author__ = 'schubert'

import itertools as itr
import warnings

import coopr.environ
from coopr.pyomo import *
from coopr.opt import SolverFactory

from Fred2.Core.Base import ACleavageSitePrediction
from Fred2.Core.Peptide import Peptide

class EpitopeAssembly(object):
    """
        Implements the epitope assembly approach proposed by Toussaint et al.
        using proteasomal cleaveage site prediction and formulationg the problem as
        TSP.

        The ILP model is implemented. So be reasonable with the size of epitope to be arranged.

        :param List(Peptide) peptides: A list of peptides which shell be arranged
        :param ACleavageSitePredictor pred: A cleavage site predictor
    """

    def __init__(self, peptides, pred, solver="glpk", verbosity=0):

        if not isinstance(pred, ACleavageSitePrediction):
            raise ValueError("Cleave site predictor must be of type ACleavageSitePrediction")

        if len(peptides) > 60:
            warnings.warn("The peptide set exceeds 60. Above this level one has to expect " +
                          "considerably long running times due to the complexity of the problem.")

        #Generate model
        #1. Generate peptides for which cleave sites have to be predicted
        #2. generate graph with dummy element
        self.__verbosity = verbosity

        peptides.append("Dummy")
        edge_matrix = {}
        fragments = {}
        seq_to_pep = {}
        for start, stop in itr.combinations(peptides, 2):
            if start == "Dummy" or stop == "Dummy":
                seq_to_pep[str(start)] = start
                seq_to_pep[str(stop)] = stop
                edge_matrix[(str(start), str(stop))] = -1
                edge_matrix[(str(stop), str(start))] = -1
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
            edge_matrix[(start, stop)] = -1.0*cleave_pred.loc[(i, len(str(start))-1), pred.name]

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

        self.instance = model.create()
        if self.__verbosity > 0:
            print "MODEL INSTANCE"
            self.instance.pprint()

    def solve(self):
        """
        Solves the Epitope Assembly problem and returns an ordered list of the peptides

        :return: list(Peptide) - An order list of the peptides (based on the string-of-beats ordering)
        """
        self.instance.preprocess()

        res = self.__solver.solve(self.instance)
        self.instance.load(res)
        if self.__verbosity > 0:
            res.write(num=1)

        return [ self.__seq_to_pep[u] for u in sorted(self.instance.u, key=lambda x: self.instance.u[x].value)]