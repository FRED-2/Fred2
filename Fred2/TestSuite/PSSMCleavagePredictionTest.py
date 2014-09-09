__author__ = 'schubert'

import unittest

from Fred2.Core.Peptide import Peptide
from Fred2.CleavagePrediction import CleavageSitePredictorFactory


class PSSMCleavagePredictonTestCase(unittest.TestCase):

    def setUp(self):
        self.seqs = [Peptide("SYFPEITHI"), Peptide("IHTIEPFYS")]

    def test_peptide_cleavage_prediction(self):
        pred = CleavageSitePredictorFactory("PCM")
        result = pred.predict(self.seqs,length=6)
        print result