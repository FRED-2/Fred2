__author__ = 'schubert'
import unittest


from Fred2.Core import Peptide, Allele
from Fred2.CleavagePrediction import CleavageSitePredictorFactory
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.EpitopeAssembly import EpitopeAssembly, ParetoEpitopeAssembly


class EpitopeAssemblyTestCase(unittest.TestCase):

    def setUp(self):
        self.peptides = [Peptide("KLLPRLPGV"), Peptide("YLYDHLAPM"), Peptide("ALYDVVSTL")]

    def test_simple_assembly(self):
        """
        Simple test if everything works. Solution manually tested for optimality.

        :return:
        """
        pred = CleavageSitePredictorFactory("PCM")
        assembler = EpitopeAssembly(self.peptides, pred, solver="cbc", verbosity=0)
        r = assembler.solve()
        self.assertEqual(r, [Peptide("KLLPRLPGV"), Peptide("ALYDVVSTL"), Peptide("YLYDHLAPM")])



    def test_pareto_assembly(self):
        cl_pred = CleavageSitePredictorFactory("PCM")
        ep_pred = EpitopePredictorFactory("SMM")
        allele = [Allele("HLA-A*02:01")]
        thresh = {a.name:10000 for a in allele}
        comp = lambda a,b: a <= b

        print ep_pred.predict(self.peptides,alleles=allele)
        #cl_pred, ep_pred, alleles, threshold, comparator, length=9

        assembler = ParetoEpitopeAssembly(self.peptides,cl_pred, ep_pred, allele, thresh, comp, solver="cbc", verbosity=1)
        r = assembler.solve(eps=1e10, order=(1,0))
        print r

        #print assembler.solve(eps=2.0)

    def test_pareto_front_assembly(self):
        cl_pred = CleavageSitePredictorFactory("PCM")
        ep_pred = EpitopePredictorFactory("SMM")
        allele = [Allele("HLA-A*02:01")]
        thresh = {a.name:10000 for a in allele}
        comp = lambda a,b: a <= b

        assembler = ParetoEpitopeAssembly(self.peptides,cl_pred, ep_pred, allele, thresh, comp, solver="cbc", verbosity=0)
        r = assembler.paretosolve()
        print r

        #print assembler.solve(eps=2.0)
