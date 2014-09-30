

"""
Unittest for ANN predictors
"""
__author__ = 'mohr, schubert'


import unittest
import pandas as pd


# Variants and Generator
import pandas
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide

#Preidctions
from Fred2.EpitopePrediction import EpitopePredictorFactory


class TestCaseANN(unittest.TestCase):

    def setUp(self):
        self.peptides = [Peptide("SYFPEITHIAAAAAA"),Peptide("IHTIEPFYSAAAAAA")]
        self.alleles = [Allele("HLA-A*02:01"),Allele("HLA-A*01:01"), Allele("HLA-A*24:02")]
        self.methods = ["NetMHCpan","NetMHC"]

    #
    # def test_netMHC_prediction(self):
    #     p = EpitopePredictorFactory("NetMHC")
    #     results = p.predict(self.peptides, self.alleles)
    #     print results
    #
    # def test_netMHCpan_prediction(self):
    #     p = EpitopePredictorFactory("NetMHCpan")
    #     results = p.predict(self.peptides, self.alleles)
    #     print results

    def test_netMHCII_prediction(self):
        p = EpitopePredictorFactory("NetMHCII")
        results = p.predict(self.peptides, [Allele("DRB1*01:01"), Allele("DRB1*04:01")])
        print results


    # def test_netMHCIIpan_prediction(self):
    #     p = EpitopePredictorFactory("NetMHCIIpan")
    #     results = p.predict(self.peptides, [Allele("DRB1*01:01"), Allele("DRB1*14:64")])
    #     print results
    #
    # def test_multiple_prediction_and_concatination(self):
    #     results = [EpitopePredictorFactory(method).predict(self.peptides, self.alleles) for method in self.methods]
    #     df1 = results[0]
    #     print df1
    #     for t in results[1:]:
    #         df1a, df2a = df1.align(t, fill_value=0)
    #         df1 = df1a+df2a
    #     print df1

if __name__ == '__main__':
    unittest.main()