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
        self.thresh = {"A*01:01":10,"B*07:02":10,"C*03:01":10}

    def test_selection_without_constraints(self):
        """
        tests if minimal selection withotu additional constraints (except the knapsack capacity) works

        #peptides obtainedn by perfroming optimization with same input and parameters by
        etk.informatik.uni-tuebingen.de/optitope

        :return:
        """
        opt = OptiTope(self.result, self.thresh, k=3, solver="cbc", verbosity=0)
        r =opt.solve()

        self.assertTrue(len(set(str(p) for p in r) - set(["GPTPLLYRL", "QYLAGLSTL", "ALYDVVSTL"])) == 0)

    def test_allele_cov_constraint(self):
        """
        tests the allele converage constraints

        :return:
        """
        #self.alleles.extend([Allele("HLA-A*02:01"),Allele("HLA-B*15:01")])
        #self.thresh.update({"A*02:01":0,"B*15:01":0})
        self.result= EpitopePredictorFactory("BIMAS").predict(self.peptides, self.alleles)
        opt = OptiTope(self.result, self.thresh, k=3, solver="cbc", verbosity=0)
        opt.activate_allele_coverage_const(0.99)
        r = opt.solve()

        self.assertTrue(len(set(str(p) for p in r) - set(["GPTPLLYRL", "QYLAGLSTL", "ALYDVVSTL"])) == 0 )

    def test_epitope_conservation_constraint(self):
        import random
        self.result = EpitopePredictorFactory("BIMAS").predict(self.peptides, self.alleles)
        conservation = {}
        for e in self.result.index.levels[0]:
            conservation[str(e)] = random.random()
        pt = OptiTope(self.result, self.thresh, k=3, solver="cbc", verbosity=0)
        pt.activate_epitope_conservation_const(0.5, conservation=conservation)
        for e in pt.solve():
            print(e, conservation[e])
