"""
Unittest for PSSM predictors
"""
__author__ = 'schubert'


import unittest


# Variants and Generator
import pandas
from Core.Allele import Allele
from Core.Peptide import Peptide

#Preidctions
from EpitopePrediction import EpitopePredictorFactory


class TestCasePSSM(unittest.TestCase):

    def setUp(self):
        self.peptides = [Peptide("SYFPEITHI"),Peptide("IHTIEPFYS")]
        self.alleles = [Allele("HLA-A*24:02"),Allele("HLA-A*02:01")]
        self.methods = ["BIMAS", "Epidemix"]

    def test_BIMAS_initialization(self):
        p = EpitopePredictorFactory("BIMAS")
        self.assertTrue(p.name == "bimas")

    def test_BIMAS_prediction(self):
        p = EpitopePredictorFactory("Epidemix")
        results = p.predict(self.peptides, self.alleles)
        print results

    def test_multiple_prediction_and_concatination(self):
        results = [EpitopePredictorFactory(method).predict(self.peptides, self.alleles) for method in self.methods[:1]]
        df1 = results[0]
        df2 = results[1]
        df1a, df2a = df1.align(df2, fill_value=0)
        df3= df1a+df2a
        print df3

if __name__ == '__main__':
    unittest.main()