__author__ = 'schubert'

import unittest

from Fred2.Core.Peptide import Peptide
from Fred2.CleavagePrediction import CleavagePredictorFactory


class PSSMCleavagePredictonTestCase(unittest.TestCase):

    def setUp(self):
        self.seqs = [Peptide("SYFPEITHI")]

    def test_peptide_cleavage_prediction(self):
        pred = CleavagePredictorFactory("PCM")
        print pred.predict(self.seqs,length=6)