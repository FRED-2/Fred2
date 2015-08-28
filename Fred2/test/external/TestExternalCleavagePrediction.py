__author__ = 'schubert'

import unittest

from Fred2.Core import Peptide
from Fred2.Core import Protein
from Fred2.Core import Transcript
from Fred2.CleavagePrediction import CleavageSitePredictorFactory


class TestExternalCleavagePredictonClass(unittest.TestCase):

    def setUp(self):
        self.seqs = [Peptide("SYFPEISYFP"),
                     Protein("IHTIEPFYSIHTIEPFYSIHTIEPFYSIHTIEPFYSIHTIEPFYS", _transcript_id="ID-01", _gene_id="FOXP3")]
        self.transcript = Transcript("")

    def test_peptide_cleavage_prediction_mixed_input(self):
        mo = CleavageSitePredictorFactory("NetChop")
        mo.predict(self.seqs)

    def test_peptide_cleavage_prediction_single_input(self):
        mo = CleavageSitePredictorFactory("NetChop")
        mo.predict(self.seqs[0])
        mo.predict(self.seqs[1])

    def test_wrong_input(self):
        with self.assertRaises(ValueError):
            mo = CleavageSitePredictorFactory("NetChop")
            mo.predict(self.seqs[0])



if __name__ == '__main__':
    unittest.main()