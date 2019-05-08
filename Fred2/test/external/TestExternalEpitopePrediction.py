"""
Unittest for external epitope prediction methods
"""

import unittest
import os

from Fred2.Core import Allele, CombinedAllele
from Fred2.Core import Peptide
from Fred2.Core import Transcript

from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.EpitopePrediction import AExternalEpitopePrediction
from Fred2.EpitopePrediction import NetMHC_3_4


#only for internal testing
class NetMHC_0_1(NetMHC_3_4):

    __version = "0.1"

    @property
    def version(self):
        return self.__version


class TestExternalEpitopePredictionClass(unittest.TestCase):

    def setUp(self):
        self.peptides_mhcI = [Peptide("SYFPEITHI"), Peptide("IHTIEPFYS")]
        self.peptides_mhcII = [Peptide("AAAAAASYFPEITHI"), Peptide("IHTIEPFYSAAAAAA")]
        self.mhcI = [Allele("HLA-B*07:02"), Allele("HLA-A*02:01")]
        self.mhcII = [Allele("HLA-DRB1*07:01"), Allele("HLA-DRB1*15:01")]
        self.mhcII_combined_alleles = [CombinedAllele("DPA1*01:03-DPB1*01:01"), CombinedAllele("DQA1*06:02-DQB1*06:31")]
        self.transcript = Transcript("")

    def test_multiple_inputs(self):
        for m in EpitopePredictorFactory.available_methods():
            for v in EpitopePredictorFactory.available_methods()[m]:
                mo = EpitopePredictorFactory(m, version=v)
                if isinstance(mo, AExternalEpitopePrediction) and not (mo.version=="0.1" and mo.name=="netmhc"):
                    print("Testing", mo.name, "version", mo.version)
                    try:
                        if any(a.name in mo.supportedAlleles for a in self.mhcII):
                            mo.predict(self.peptides_mhcII, alleles=self.mhcII)
                        if any(a.name in mo.supportedAlleles for a in self.mhcII_combined_alleles):
                            mo.predict(self.peptides_mhcII, alleles=self.mhcII_combined_alleles)
                        if any(a.name in mo.supportedAlleles for a in self.mhcI):
                            mo.predict(self.peptides_mhcI, alleles=self.mhcI)
                        print("Success")
                    except RuntimeError as e: #catch only those stemming from binary unavailability
                        if "could not be found in PATH" not in e.message:
                            raise e #all others do not except
                        else:
                            print(mo.name, "not available")

    def test_single_epitope_input(self):
        for m in EpitopePredictorFactory.available_methods():
            for v in EpitopePredictorFactory.available_methods()[m]:
                mo = EpitopePredictorFactory(m, version=v)
                if isinstance(mo, AExternalEpitopePrediction) and not (mo.version=="0.1" and mo.name=="netmhc"):
                    print("Testing", mo.name, "version", mo.version)
                    try:
                        if any(a.name in mo.supportedAlleles for a in self.mhcII):
                            mo.predict(self.peptides_mhcII[0], alleles=self.mhcII)
                        if any(a.name in mo.supportedAlleles for a in self.mhcII_combined_alleles):
                            mo.predict(self.peptides_mhcII[0], alleles=self.mhcII_combined_alleles)
                        if any(a.name in mo.supportedAlleles for a in self.mhcI):
                            mo.predict(self.peptides_mhcI[0], alleles=self.mhcI)
                        print("Success")
                    except RuntimeError as e: #catch only those stemming from binary unavailability
                        if "could not be found in PATH" not in e.message:
                            raise e #all others do not except
                        else:
                            print(mo.name, "not available")

    def test_single_allele_input(self):
        for m in EpitopePredictorFactory.available_methods():
            for v in EpitopePredictorFactory.available_methods()[m]:
                mo = EpitopePredictorFactory(m, version=v)
                if isinstance(mo, AExternalEpitopePrediction) and not (mo.version=="0.1" and mo.name=="netmhc"):
                    print("Testing", mo.name, "version", mo.version)
                    try:
                        if any(a.name in mo.supportedAlleles for a in self.mhcII):
                            mo.predict(self.peptides_mhcII, alleles=self.mhcII[0])
                        if any(a.name in mo.supportedAlleles for a in self.mhcII_combined_alleles):
                            mo.predict(self.peptides_mhcII, alleles=self.mhcII_combined_alleles[0])
                        if any(a.name in mo.supportedAlleles for a in self.mhcI):
                            mo.predict(self.peptides_mhcI, alleles=self.mhcI[0])
                        print("Success")
                    except RuntimeError as e: #catch only those stemming from binary unavailability
                        if "could not be found in PATH" not in e.message:
                            raise e #all others do not except
                        else:
                            print(mo.name, "not available")

    def test_wrong_epitope_input(self):
        with self.assertRaises(ValueError):
            EpitopePredictorFactory("NetMHC").predict(self.transcript, alleles=self.mhcI)

    def test_wrong_allele_input(self):
        with self.assertRaises(ValueError):
            EpitopePredictorFactory("NetMHC").predict(self.mhcI, alleles=self.transcript)

    def test_wrong_internal_to_external_version(self):
        with self.assertRaises(RuntimeError):
            EpitopePredictorFactory("NetMHC", version="0.1").predict(self.peptides_mhcI, alleles=self.mhcI)

    def test_path_option_and_optional_parameters_netmhc(self):
        netmhc = EpitopePredictorFactory("NetMHC")
        exe = netmhc.command.split()[0]
        for try_path in os.environ["PATH"].split(os.pathsep):
            try_path = try_path.strip('"')
            exe_try = os.path.join(try_path, exe).strip()
            if os.path.isfile(exe_try) and os.access(exe_try, os.X_OK):
                r = netmhc.predict(self.peptides_mhcI, alleles=self.mhcI, command=exe_try, options="--sort", chunksize=1)
                self.assertTrue(len(r) == len(self.peptides_mhcI))
                self.assertAlmostEqual(r["A*02:01"]["SYFPEITHI"]["netmhc"], 0.150579105869, places=7, msg=None, delta=None)
                self.assertAlmostEqual(r["A*02:01"]["IHTIEPFYS"]["netmhc"], 0.0619540879359, places=7, msg=None, delta=None)

    def test_path_and_optional_parameters_netctl(self):
        netctlpan = EpitopePredictorFactory("NetCTLpan")
        exe = netctlpan.command.split()[0]
        for try_path in os.environ["PATH"].split(os.pathsep):
            try_path = try_path.strip('"')
            exe_try = os.path.join(try_path, exe).strip()
            if os.path.isfile(exe_try) and os.access(exe_try, os.X_OK):
                print(netctlpan.predict(self.peptides_mhcI, alleles=self.mhcI,
                                        commad=exe_try,
                                        options="-wt 0.05 -wc 0.225 -ethr 0.5"))

if __name__ == '__main__':
    unittest.main()