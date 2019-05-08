"""
Unittest of versioning system
"""

import unittest

# Variants and Generator
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide

#Preidctions
from Fred2.EpitopePrediction import EpitopePredictorFactory,BIMAS, AEpitopePrediction
from Fred2.TAPPrediction import TAPPredictorFactory
from Fred2.CleavagePrediction import CleavageFragmentPredictorFactory, CleavageSitePredictorFactory


class BIMAS2(BIMAS):
    __version = "2.0"

    @property
    def version(self):
        return self.__version


class TestCaseEpitopePrediction(unittest.TestCase):

    def setUp(self):
        #Peptides of different length 9,10,11,12,13,14,15
        self.peptides_mhcI = [Peptide("SYFPEITHI"), Peptide("IHTIEPFYS")]
        self.peptides_fragment = [Peptide("IHTIEPFYSAA")]
        self.mhcI = [Allele("HLA-B*15:01"), Allele("HLA-A*02:01")]
        self.mhcII = [Allele("HLA-DRB1*07:01"), Allele("HLA-DRB1*15:01")]

    def test_epitope_prediction_specific_version(self):
        print(EpitopePredictorFactory("BIMAS", version="1.0").predict(self.peptides_mhcI, self.mhcI))

    def test_epitope_prediction_no_version(self):
        print(EpitopePredictorFactory("BIMAS").predict(self.peptides_mhcI, self.mhcI))

    def test_epitope_prediction_available_methods(self):
        print(EpitopePredictorFactory.available_methods())

    def test_multiple_predictors_names_different_version(self):
        self.assertTrue(EpitopePredictorFactory("BIMAS", version="1.0").version == "1.0")
        self.assertTrue(EpitopePredictorFactory("BIMAS", version="2.0").version == "2.0")

    @unittest.expectedFailure
    def test_epitope_prediction_unsupported_version(self):
        print(EpitopePredictorFactory("BIMAS", version="4.0").predict(self.peptides_mhcI, self.mhcI))

    def test_TAP_prediction_specific_version(self):
        print(TAPPredictorFactory("SVMTAP", version="1.0").predict(self.peptides_mhcI))

    def test_TAP_prediction_no_version(self):
        print(TAPPredictorFactory("SVMTAP").predict(self.peptides_mhcI))

    def test_TAP_prediction_available_methods(self):
        print(TAPPredictorFactory.available_methods())

    @unittest.expectedFailure
    def test_TAP_prediction_unsupported_version(self):
        print(TAPPredictorFactory("SVMTAP", version="5.0").predict(self.peptides_mhcI))

    def test_CleavageSite_prediction_specific_version(self):
        print(CleavageSitePredictorFactory("PCM", version="1.0").predict(self.peptides_mhcI))

    def test_CleavageSite_prediction_no_version(self):
        print(CleavageSitePredictorFactory("PCM").predict(self.peptides_mhcI))

    def test_CleavageSite_prediction_available_methods(self):
        print(CleavageSitePredictorFactory.available_methods())

    @unittest.expectedFailure
    def test_CleavageSite_prediction_unsupported_version(self):
        print(CleavageSitePredictorFactory("PCM", version="2341.0").predict(self.peptides_mhcI))

    def test_CleavageFrag_prediction_specific_version(self):
        print(CleavageFragmentPredictorFactory("Ginodi", version="1.0").predict(self.peptides_fragment))

    def test_CleavageFrag_prediction_no_version(self):
        print(CleavageFragmentPredictorFactory("Ginodi").predict(self.peptides_fragment))

    def test_CleavageFrag_prediction_available_methods(self):
        print(CleavageFragmentPredictorFactory.available_methods())

    @unittest.expectedFailure
    def test_CleavageFrag_prediction_unsupported_version(self):
        print(CleavageFragmentPredictorFactory("Ginodi", version="1.234").predict(self.peptides_fragment))

if __name__ == '__main__':
    unittest.main()