import unittest


from Fred2.Core.Peptide import Peptide
from Fred2.TAPPrediction import TAPePredictorFactory

class TAPPredictionTestCaste(unittest.TestCase):

    def setUp(self):
        self.peptides = [Peptide("SYFPEITHI"),Peptide("IHTIEPFYS")]


    def test_svmtap_simple(self):
        """
            Tests SVMTAP prediction

            not compared yet (dont know where)
        """

        pred = TAPePredictorFactory("SVMTAP")
        r = pred.predict(self.peptides)
        print r