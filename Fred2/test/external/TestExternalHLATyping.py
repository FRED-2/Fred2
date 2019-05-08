"""
Unittest for external HLA typing  methods
"""
import sys
sys.path.append("/home/schubert/Dropbox/PhD/Porgramming/Fred2/")
import unittest
import os
import subprocess
import pprint
from Fred2.Core import Allele
from Fred2.HLAtyping import HLATypingFactory

from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.CleavagePrediction import CleavageSitePredictorFactory, CleavageFragmentPredictorFactory
from Fred2.TAPPrediction import TAPPredictorFactory

class TestExternalHLATypingClass(unittest.TestCase):

    def setUp(self):
        self.file1 = "/Users/schubert/Dropbox/PhD/Porgramming/OptiType/test/rna/CRC_81_N_1_fished.fastq"
        self.file2 = "/Users/schubert/Dropbox/PhD/Porgramming/OptiType/test/rna/CRC_81_N_2_fished.fastq"

    def test_optitype(self):
        os.environ["PATH"] += os.pathsep + "/Users/schubert/Dropbox/PhD/Porgramming/OptiType/"
        #print os.environ["PATH"]
        opti = HLATypingFactory("OptiType")
        print(opti.predict("/Users/schubert/Dropbox/PhD/Porgramming/OptiType/test/exome/NA11995_SRR766010_1_fished.fastq", "/tmp/", options="-d"))

    def test_seq2HLA(self):
        origin = "/home/schubert/Dropbox/PhD/software/seq2hla"
        os.environ["PATH"] += os.pathsep + origin
        seq2HLA = HLATypingFactory("Seq2HLA")
        print(seq2HLA.predict("/home/schubert/Desktop/ERR009105_1.fastq.gz", origin+"/delete", options="-2 /home/schubert/Desktop/ERR009105_2.fastq.gz"))

    def test_athlates(self):
        origin = "/Users/schubert/Dropbox/PhD/software/Athlates_2014_04_26/bin"
        os.environ["PATH"] += os.pathsep + origin
        atlates = HLATypingFactory("athlates")
        print(atlates.predict("/Users/schubert/Dropbox/PhD/software/Athlates_2014_04_26/demo/HG01756/HG01756_a.sort.bam",
                              "/Users/schubert/Dropbox/PhD/software/Athlates_2014_04_26/demo/output/Fred2_",
                              options="-msa /Users/schubert/Dropbox/PhD/software/Athlates_2014_04_26/db/msa/A_nuc.txt",
                              delete=False))


    def test_polysolver(self):
        command = ['source /Users/schubert/Dropbox/phd/software/polysolver/scripts/config.bash && env']
        proc = subprocess.Popen(command, stdout = subprocess.PIPE,shell=True)
        for line in proc.stdout:
            (key, _, value) = line.partition("=")
            os.environ[key] = value.strip()
        proc.communicate()
        origin = "/Users/schubert/Dropbox/phd/software/polysolver/scripts/"
        os.environ["PATH"] += os.pathsep + origin
        pprint.pprint(dict(os.environ))



        command = ["/Users/schubert/Dropbox/phd/software/polysolver/debugg.sh"]
        proc = subprocess.Popen(command, stdout = subprocess.PIPE,shell=True)

        print("SAMTool call")
        print("\n".join( l for l in proc.stdout))
        print()
        print(str(proc.stderr))
        proc.communicate()
        #sys.exit()
        polysolver = HLATypingFactory("polysolver")
        print(polysolver.predict("/Users/schubert/Dropbox/phd/software/polysolver/test/test.bam",
                                 "/Users/schubert/Dropbox/phd/software/polysolver/debugging2/", options="Unknown 1 hg19 STDFQ 0",
                                 delete=False,
                                 command=origin+"shell_call_hla_type"))
    #


if __name__ == '__main__':
    unittest.main()

