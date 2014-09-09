__author__ = 'Schubert'

import unittest


from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.EpitopeSelection.OptiTope import OptiTope

class OptiTopeTestCase(unittest.TestCase):
    """
        Unittest for OptiTope
    """

    def setUp(self):
        self.proteins=[]
        self.alleles = [Allele("HLA-A*01:01"),Allele("HLA-B*07:02"), Allele("HLA-C*03:01")]
        self.peptides = [Peptide(p) for p in """SFSIFLLAL
GHRMAWDMM
VYEADDVIL
CFTPSPVVV
FLLLADARV
GPADGMVSK
YLYDHLAPM
GLRDLAVAV
GPTPLLYRL
TWVLVGGVL
IELGGKPAL
LAGGVLAAV
QYLAGLSTL
NFVSGIQYL
VLSDFKTWL
ARPDYNPPL
KLLPRLPGV
RHTPVNSWL
GLYLFNWAV
ALYDVVSTL
RRCRASGVL
WPLLLLLLA
VTYSLTGLW
YFVIFFVAA""".split()]
        self.result= EpitopePredictorFactory("BIMAS").predict(self.peptides, self.alleles)
        self.thresh = {"A*01:01":0,"B*07:02":0,"C*03:01":0}
    #
    # def test_selection_without_constraints(self):
    #     """
    #     tests if minimal selection withotu additional constraints (except the knapsack capacity) works
    #
    #     #peptides obtainedn by perfroming optimization with same input and parameters by
    #     etk.informatik.uni-tuebingen.de/optitope
    #
    #     :return:
    #     """
    #     opt = OptiTope(self.result, self.alleles,self.thresh,k=3,solver="cplex",verbosity=0)
    #     r =opt.solve()
    #     self.assertTrue(len(set(str(p) for p in r) - set(["GPTPLLYRL", "QYLAGLSTL", "ALYDVVSTL"])) == 0 )

    def test_allele_cov_constraint(self):
        """
        tests the allele converage constraints

        :return:
        """
        self.alleles.extend([Allele("HLA-A*02:01"),Allele("HLA-B*15:01")])
        self.thresh.update({"A*02:01":0,"B*15:01":0})
        self.result= EpitopePredictorFactory("BIMAS").predict(self.peptides, self.alleles)
        opt = OptiTope(self.result, self.alleles,self.thresh,k=3,solver="cplex",verbosity=1)
        opt.activate_allele_coverage_const(0.8)
        r = opt.solve()
        res_df = self.result.xs(self.result.index.values[0][1], level="Method")
        peps = [str(p) for p in r]

        probs = {"A*01:01":0.5, "A*02:01":0.5, "B*07:02":0.5,"B*15:01":0.5,"C*03:01":1.0}
        res_df = res_df.loc[peps, :]
        res_df = res_df[[a.name for a in self.alleles]]
        res_df = res_df[res_df.apply(lambda x: any(x[a.name] > self.thresh[a.name] for a in self.alleles), axis=1)]

        print res_df.apply(lambda x: sum( x[c]*probs[c] for c in res_df.columns),axis=1)

        self.assertTrue(len(set(str(p) for p in r) - set(["ALYDVVSTL", "KLLPRLPGV", "YLYDHLAPM"])) == 0 )

