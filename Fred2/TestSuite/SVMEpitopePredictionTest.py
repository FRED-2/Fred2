import unittest
import pandas as pd

from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory

class SVMEpitopePredictionTestCase(unittest.TestCase):

    def setUp(self):
        self.peptides = [Peptide("SYFPEITHI"),Peptide("IHTIEPFYS")]
        self.alleles = [Allele("HLA-A*24:02"),Allele("HLA-A*02:01")]

    def test_svmhc_simple(self):
        """
        Output compared to SVMHC server output.
        """
        pred = EpitopePredictorFactory("SVMHC")
        r = pred.predict(self.peptides, self.alleles)
        print r

    def test_unitope_simple(self):
        """
        Classification in concordance with webservice
        """
        pred = EpitopePredictorFactory("UniTope")
        r = pred.predict(self.peptides, self.alleles)
        print r

    def test_MHCIIMulti_simple(self):
        pred = EpitopePredictorFactory("MHCIIMulti")
        r = pred.predict(self.peptides, [Allele("DRB1*07:09")])
        print r