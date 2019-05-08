__author__ = 'schubert'

import unittest
import os

from Fred2.Core import Peptide
from Fred2.Core import Protein
from Fred2.Core import Transcript
from Fred2.CleavagePrediction import CleavageSitePredictorFactory


class TestExternalCleavagePredictonClass(unittest.TestCase):

    def setUp(self):
        self.seqs = [Peptide("SYFPEISYFP"),
                     Protein("IHTIEPFYSIHTIEPFYSIHTIEPFYSIHTIEPFYSIHTIEPFYS", transcript_id="ID-01", gene_id="FOXP3")]
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

    def test_path_option_and_optionl_parameters(self):
        netchop = CleavageSitePredictorFactory("NetChop")
        exe = netchop.command.split()[0]
        for try_path in os.environ["PATH"].split(os.pathsep):
            try_path = try_path.strip('"')
            exe_try = os.path.join(try_path, exe).strip()
            if os.path.isfile(exe_try) and os.access(exe_try, os.X_OK):
                print(exe_try)
                netchop.predict(self.seqs, path=exe_try, options="-v 1")

if __name__ == '__main__':
    unittest.main()