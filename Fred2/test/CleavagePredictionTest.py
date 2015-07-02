__author__ = 'schubert'

import unittest

from Fred2.Core.Peptide import Peptide
from Fred2.Core.Protein import Protein
from Fred2.CleavagePrediction import CleavageSitePredictorFactory, CleavageFragmentPredictorFactory


class PSSMCleavagePredictonTestCase(unittest.TestCase):

    def setUp(self):
        self.seqs = [Peptide("SYFPEISYFP"), Protein("IHTIEPFYSIHTIEPFYSIHTIEPFYSIHTIEPFYSIHTIEPFYS","ID-01","FOXP3")]
        self.fragments= [Peptide("FSYFPEITHIR"), Peptide("FIHTIEPFYSR")]

    def test_peptide_cleavage_prediction_mixed_input(self):
        for m in CleavageSitePredictorFactory.available_methods():
            if m != "NetChop":
                mo = CleavageSitePredictorFactory(m)
                mo.predict(self.seqs)

    def test_peptide_cleavage_prediction_single_input(self):
        for m in CleavageSitePredictorFactory.available_methods():
            if m != "NetChop":
                mo = CleavageSitePredictorFactory(m)
                mo.predict(self.seqs[0])
                mo.predict(self.seqs[1])

    def test_cleavage_fragment_prediction_multiple_input(self):
        for m in CleavageFragmentPredictorFactory.available_methods():
            pred = CleavageFragmentPredictorFactory(m)
            pred.predict(self.fragments)

    def test_cleavage_fragment_prediction_single_input(self):
        for m in CleavageFragmentPredictorFactory.available_methods():
            pred = CleavageFragmentPredictorFactory(m)
            pred.predict(self.fragments[0])


if __name__ == '__main__':
    unittest.main()