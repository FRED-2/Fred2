"""
Unittest for PSSM predictors
"""
__author__ = 'schubert'


import unittest

# Variants and Generator
from Fred2.Core import Allele
from Fred2.Core import Peptide

#Preidctions
from Fred2.EpitopePrediction import EpitopePredictorFactory, AExternalEpitopePrediction


class TestCaseEpitopePrediction(unittest.TestCase):

    def setUp(self):
        #Peptides of different length 9,10,11,12,13,14,15
        self.peptides_mhcI = [Peptide("SYFPEITHI"), Peptide("IHTIEPFYS")]
        self.peptides_mhcII = [Peptide("SYFPEITHI"), Peptide("IHTIEPFYSAAAAAA")]
        self.mhcI = [Allele("HLA-B*15:01"), Allele("HLA-A*02:01")]
        self.mhcII = [Allele("HLA-DRB1*07:01"), Allele("HLA-DRB1*15:01")]

    def test_multiple_peptide_input_mhcI(self):
            for m in EpitopePredictorFactory.available_methods():
                model = EpitopePredictorFactory(m)
                if not isinstance(model, AExternalEpitopePrediction):
                    if all(a.name in model.supportedAlleles for a in self.mhcI):
                        res = model.predict(self.peptides_mhcI, alleles=self.mhcI)

    def test_single_peptide_input_mhcI(self):
            for m in EpitopePredictorFactory.available_methods():
                model = EpitopePredictorFactory(m)
                if not isinstance(model, AExternalEpitopePrediction):
                    if all(a.name in model.supportedAlleles for a in self.mhcI):
                        res = model.predict(self.peptides_mhcI[0], alleles=self.mhcI[0])

    def test_multiple_peptide_input_mhcII(self):
            for m in EpitopePredictorFactory.available_methods():
                model = EpitopePredictorFactory(m)
                if not isinstance(model, AExternalEpitopePrediction):
                    if all(a.name in model.supportedAlleles for a in self.mhcII) and m != "MHCIIMulti":
                        res = model.predict(self.peptides_mhcII, alleles=self.mhcII)

    def test_single_peptide_input_mhcII(self):
            for m in EpitopePredictorFactory.available_methods():
                model = EpitopePredictorFactory(m)
                if not isinstance(model, AExternalEpitopePrediction):
                    if all(a.name in model.supportedAlleles for a in self.mhcII):
                        res = model.predict(self.peptides_mhcII[0], alleles=self.mhcII[1])


if __name__ == '__main__':
    unittest.main()