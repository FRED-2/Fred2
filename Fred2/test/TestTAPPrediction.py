import unittest


from Fred2.Core.Peptide import Peptide
from Fred2.TAPPrediction import TAPPredictorFactory

class TAPPredictionTestCaste(unittest.TestCase):

    def setUp(self):
        self.peptides = [Peptide("SYFPEITHI"), Peptide("IHTIEPFYS")]

    def test_tap_multiple_peptide_input(self):
        """
            Tests SVMTAP prediction

            not compared yet (dont know where)
        """
        for m in TAPPredictorFactory.available_methods():
            pred = TAPPredictorFactory(m)
            r = pred.predict(self.peptides)
            print(r)

    def test_tap_single_peptide_input(self):
        """
            Tests SVMTAP prediction

            not compared yet (dont know where)
        """
        for m in TAPPredictorFactory.available_methods():
            pred = TAPPredictorFactory(m)
            r = pred.predict(self.peptides[0])
            print(r)

    def test_smmtap_abitrary_peptide_length(self):
        smmtap = TAPPredictorFactory("smmtap")
        peptides = [Peptide("SYFPEITHI"), Peptide("IHTIEPFYSA"), Peptide("IHTIEPFYSAA")]
        print(smmtap.predict(peptides))

    def test_peptide_chunksize(self):
        for m in TAPPredictorFactory.available_methods():
            pred = TAPPredictorFactory(m)
            r = pred.predict(self.peptides, chunksize=1)
            print(r)