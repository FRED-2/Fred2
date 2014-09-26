"""
Unittest for PSSM predictors
"""
__author__ = 'schubert'


import unittest
import pandas as pd
import numpy as np

# Variants and Generator
import pandas
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide

#Preidctions
from Fred2.EpitopePrediction import EpitopePredictorFactory


class TestCasePSSM(unittest.TestCase):

    def setUp(self):
        self.peptides = [Peptide("SYFPEITHI"),Peptide("IHTIEPFYS")]
        self.alleles = [Allele("HLA-A*24:02"),Allele("HLA-A*02:01"),Allele("HLA-DRB1*01:01")]
        self.methods = ["BIMAS","Epidemix","Syfpeithi", "SMM", "SMMPMBEC", "ARB", "ComblibSidney2008","TEPITOPEpan"]

    def test_BIMAS_initialization(self):
        p = EpitopePredictorFactory("BIMAS")
        self.assertTrue(p.name == "bimas")

    def test_BIMAS_prediction(self):
        p = EpitopePredictorFactory("Epidemix")
        results = p.predict(self.peptides, self.alleles)
        print results

    def test_multiple_prediction_and_concatination(self):
        results = [EpitopePredictorFactory(method).predict(self.peptides, self.alleles) for method in self.methods]
        df1 = results[0]

        for r in results[1:]:
            tmp_a,tmp_b = df1.align(r, fill_value=0)
            df1 = tmp_a+tmp_b

        print df1
        #for t in df1:
        #    print type(t)

        #df=df3.xs(df3.index.values[0][1], level="Method")
        #for p in df.itertuples():
        #    for a,s in zip(df.columns,p[1:]):
        #        print p[0], a, s
        #        print type(p[0])

        #print df3

if __name__ == '__main__':
    unittest.main()