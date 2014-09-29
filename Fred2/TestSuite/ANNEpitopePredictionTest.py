

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
        self.peptides = [Peptide("SYFPEITHI"),Peptide("IHTIEPFYS")]
        self.alleles = [Allele("HLA-A*02:01"),Allele("HLA-A*01:01"), Allele("HLA-A*23:18")]
        self.methods = ["NetMHC"]

    def test_netMHC_initialization(self):
        p = EpitopePredictorFactory('NetMHCpan')
        r = p.predict(self.peptides, self.alleles)
        print r
        #self.assertTrue(p.name == 'netmhc')

    # def test_netMHC_prediction(self):
    #     p = EpitopePredictorFactory("NetMHC")
    #     results = p.predict(self.peptides, self.alleles)
    #     print results
    #
    # def test_netMHCpan_initialization(self):
    #     p = EpitopePredictorFactory("NetMHCpan")
    #     self.assertTrue(p.name == "netmhcpan")
    #
    # def test_netMHC_prediction(self):
    #     p = EpitopePredictorFactory("NetMHCpan")
    #     results = p.predict(self.peptides, self.alleles)
    #     print results
    #
    # def test_multiple_prediction_and_concatination(self):
    #     results = [EpitopePredictorFactory(method).predict(self.peptides, self.alleles) for method in self.methods]
    #     df1 = results[0]
    #     print df1
    #     #for t in df1:
    #     #    print type(t)
    #     #df2 = results[1]
    #     #df1a, df2a = df1.align(df2, fill_value=0)
    #     #df3= df1a+df2a
    #     #print df3

        #print df3
if __name__ == '__main__':
    unittest.main()