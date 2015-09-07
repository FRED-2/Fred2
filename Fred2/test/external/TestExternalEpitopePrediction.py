"""
Unittest for external epitope prediction methods
"""



import unittest

from Fred2.Core import Allele
from Fred2.Core import Peptide
from Fred2.Core import Transcript

from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.EpitopePrediction import AExternalEpitopePrediction


class TestExternalEpitopePredictionClass(unittest.TestCase):

    def setUp(self):
        self.peptides_mhcI = [Peptide("SYFPEITHI"),Peptide("IHTIEPFYS")]
        self.peptides_mhcII = [Peptide("AAAAAASYFPEITHI"),Peptide("IHTIEPFYSAAAAAA")]
        self.mhcI = [Allele("HLA-B*15:01"),Allele("HLA-A*02:01")]
        self.mhcII = [Allele("HLA-DRB1*07:01"), Allele("HLA-DRB1*15:01")]
        self.transcript = Transcript("")

    def test_multiple_inputs(self):
        for m in EpitopePredictorFactory.available_methods():
            mo = EpitopePredictorFactory(m)
            if isinstance(mo, AExternalEpitopePrediction):
                if any(a.name in mo.supportedAlleles for a in self.mhcII):
                    mo.predict(self.peptides_mhcII,alleles=self.mhcII)
                else:
                    mo.predict(self.peptides_mhcI,alleles=self.mhcI)

    def test_single_epitope_input(self):
        for m in EpitopePredictorFactory.available_methods():
            mo = EpitopePredictorFactory(m)
            if isinstance(mo, AExternalEpitopePrediction):
                if any(a.name in mo.supportedAlleles for a in self.mhcII):
                    mo.predict(self.peptides_mhcII[0],alleles=self.mhcII)
                else:
                    mo.predict(self.peptides_mhcI[0],alleles=self.mhcI)

    def test_single_allele_input(self):
        for m in EpitopePredictorFactory.available_methods():
            mo = EpitopePredictorFactory(m)
            if isinstance(mo, AExternalEpitopePrediction):
                if any(a.name in mo.supportedAlleles for a in self.mhcII):
                    mo.predict(self.peptides_mhcII, alleles=self.mhcII[0])
                else:
                    mo.predict(self.peptides_mhcI, alleles=self.mhcI[0])

    def test_wrong_epitope_input(self):
        with self.assertRaises(ValueError):
            EpitopePredictorFactory("NetMHC").predict(self.transcript, alleles=self.mhcI)

    def test_wrong_allele_input(self):
        with self.assertRaises(ValueError):
            EpitopePredictorFactory("NetMHC").predict(self.mhcI, alleles=self.transcript)