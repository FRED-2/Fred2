"""
Unittest for PSSM predictors
"""
__author__ = 'schubert'


import unittest


# Variants and Generator
from Core.Allele import Allele
from Core.Peptide import Peptide

#Preidctions
from Prediction import EpitopePredictorFactory


class TestCasePSSM(unittest.TestCase):

    def setUp(self):
        self.peptides = [Peptide("SYFPEITHI")]
        self.alleles = [Allele("HLA-A*31:01")]

    def test_BIMAS_initialization(self):
        p = EpitopePredictorFactory("BIMAS")
        self.assertTrue(p.name == "bimas")

    def test_BIMAS_prediction(self):
        p = EpitopePredictorFactory("BIMAS")
        results = p.predict(self.peptides, self.alleles)
        print results



if __name__ == '__main__':
    unittest.main()