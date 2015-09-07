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
from Fred2.EpitopePrediction import AExternalEpitopePrediction


class TestCaseEpitopePrediction(unittest.TestCase):

    def setUp(self):
        #Peptides of different length 9,10,11,12,13,14,15
        self.peptides_mhcI = [Peptide("SYFPEITHI"), Peptide("IHTIEPFAA")]
        self.peptides_mhcII = [Peptide("SYFPEITHI"), Peptide("IHTIEPFYSAAAAAA")]
        self.mhcI = [Allele("HLA-B*15:01"), Allele("HLA-A*02:01")]
        self.mhcII = [Allele("HLA-DRB1*07:01"), Allele("HLA-DRB1*15:01")]


    def test_multiple_peptide_input_mhcI(self):
            for m in EpitopePredictorFactory.available_methods():
                model = EpitopePredictorFactory(m)
                if not isinstance(model, AExternalEpitopePrediction):
                    if any(a.name in model.supportedAlleles for a in self.mhcI):
                        res = model.predict(self.peptides_mhcI, alleles=self.mhcI)

    def test_single_peptide_input_mhcI(self):
            for m in EpitopePredictorFactory.available_methods():
                model = EpitopePredictorFactory(m)
                if not isinstance(model, AExternalEpitopePrediction):
                    if any(a.name in model.supportedAlleles for a in self.mhcI):
                        res = model.predict(self.peptides_mhcI[0], alleles=self.mhcI[1])

    def test_multiple_peptide_input_mhcII(self):
            for m in EpitopePredictorFactory.available_methods():
                model = EpitopePredictorFactory(m)
                if not isinstance(model, AExternalEpitopePrediction):
                    if any(a.name in model.supportedAlleles for a in self.mhcII) and m != "MHCIIMulti":
                        res = model.predict(self.peptides_mhcII, alleles=self.mhcII)


    def test_single_peptide_input_mhcII(self):
            for m in EpitopePredictorFactory.available_methods():
                model = EpitopePredictorFactory(m)
                if not isinstance(model, AExternalEpitopePrediction):
                    if any(a.name in model.supportedAlleles for a in self.mhcII[:1]):
                        print model.name
                        res = model.predict(self.peptides_mhcII[0], alleles=self.mhcII[0])

if __name__ == '__main__':
    unittest.main()