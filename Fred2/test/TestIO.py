from unittest import TestCase
import copy

from Fred2.Core import Peptide
from Fred2.Core import Allele
from Fred2.IO import FileReader
from Fred2.Core import Transcript
from Fred2.Core import Variant
from Fred2.Core import VariationType
from Fred2.Core import MutationSyntax

__author__ = 'walzer'


class TestIO(TestCase):

    def setUp(self):
        self.ale_path = "../Data/examples/alleles.txt"
        self.ale_zonk_path = "../Data/examples/alleles_defect.txt"
        self.ale_no_path = "../Data/examples/magic.txt"
        self.fa_path = "../Data/examples/testSequences.fasta"
        self.fa_unconventional_path = "../Data/examples/testSequences_d2s.fasta"
        self.ano_path = "../tutorials/data/test_annovar.out"

    def test_read_lines(self):
        alleles = FileReader.read_lines(self.ale_path, in_type=Allele)
        self.assertEqual(len(alleles), 2)
        self.assertRaises(IOError, FileReader.read_lines, self.ale_no_path, in_type=Allele)
        self.assertRaises(ValueError, FileReader.read_lines, self.ale_zonk_path, in_type=Allele)

    def test_read_fasta(self):
        seqs = FileReader.read_fasta(self.fa_path)
        self.assertEqual(len(seqs), 2)
        self.assertRaises(IndexError, FileReader.read_fasta, self.fa_unconventional_path) # no "|"

    def test_read_annovar_exonic(self):
        ano = FileReader.read_annovar_exonic(self.ano_path)
        self.assertEqual(len(ano), 5)
