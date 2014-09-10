__author__ = 'schubert'
import unittest


from Fred2.Core.Peptide import Peptide
from Fred2.CleavagePrediction import CleavageSitePredictorFactory
from Fred2.EpitopeAssembly.EpitopeAssembly import EpitopeAssembly


class EpitopeAssemblyTestCase(unittest.TestCase):

    def setUp(self):
        self.peptides = [Peptide("KLLPRLPGV"), Peptide("YLYDHLAPM"), Peptide("ALYDVVSTL")]


    def test_simple_assembly(self):
        pred = CleavageSitePredictorFactory("PCM")
        assembler = EpitopeAssembly(self.peptides,pred,solver="cplex",verbosity=1)
        print assembler.solve()