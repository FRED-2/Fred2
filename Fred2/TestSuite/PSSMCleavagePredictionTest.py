__author__ = 'schubert'

import unittest

from Fred2.Core.Peptide import Peptide
from Fred2.CleavagePrediction import CleavageSitePredictorFactory, CleavageFragmentPredictorFactory


class PSSMCleavagePredictonTestCase(unittest.TestCase):

    def setUp(self):
        self.seqs = [Peptide("SYFPEI"), Peptide("IHTIEPFYS")]
        self.fragments= [Peptide("FSYFPEITHIR"), Peptide("FIHTIEPFYSR")]

    def test_peptide_cleavage_prediction(self):
        pred = CleavageSitePredictorFactory("PCM")
        result = pred.predict(self.seqs, length=6)
        print result
        rs= result.xs((slice(None),pred.cleavagePos-1))
        print "score ", rs.loc["seq_0","pcm"]
        s = result.ix["seq_0"]["Seq"]
        print "astype()",s.astype(str)
        p = "".join(s)
        print type(p),p
        print s.index

    def test_cleavage_fragment_prediction(self):
        pred = CleavageFragmentPredictorFactory("PSSMGinodi")
        result = pred.predict(self.fragments)
        print result