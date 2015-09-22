from unittest import TestCase
import copy

from Fred2.Core import Peptide
from Fred2.Core import Protein
from Fred2.Core import Transcript
from Fred2.Core import Variant
from Fred2.Core import VariationType
from Fred2.Core import MutationSyntax

__author__ = 'walzer'


class TestPeptide(TestCase):

    def setUp(self):
        self.simple = Peptide("SYFPEITHI")

        self.gcg_ps = "MKSIYFVAGLFVMLVQGSWQRSLQDTEEKSRSFSASQADPLSDPDQMNEDKRHSQGTFTSDYSKYLDSRRAQDFVQWLMNTKRNRNNIAKRHDEFERHAEGTFTSDVSSYLEGQAAKEFIAWLVKGRGRRDFPEEVAIVEELGRRHADGSFSDEMNTILDNLAARDFINWLIQTKITDRK"
        self.gcg_t1 = Transcript("", _transcript_id="GLUC_HUMAN")
        gcg_p1 = Protein(self.gcg_ps, _transcript_id='GLUC_HUMAN', _orig_transcript=self.gcg_t1)
        self.w_p = Peptide("PROTEIN", {gcg_p1:[0]})
        self.gcg_p1 = gcg_p1
        self.gcg_v1 = Variant("rs5650", VariationType.SNP, 2, 162145588, 'G', 'T',
                         {"GLUC_HUMAN": MutationSyntax("GLUC_HUMAN", 344, 115, "c.344C>A", "p.A115D")}, False, False)
        gcg_p1_copy = copy.deepcopy(gcg_p1)
        gcg_p1_copy.vars = {0:[self.gcg_v1]}
        self.w_v = Peptide("VARIANT", {gcg_p1_copy:[0]})

        #self.gcg_ts = "gcatagaatgcagatgagcaaagtgagtgggagagggaagtcatttgtaacaaaaactcattatttacagatgagaaatttatattgtcagcgtaatatctgtgaggctaaacagagctggagagtatataaaagcagtgcgccttggtgcagaagtacagagcttaggacacagagcacatcaaaagttcccaaagagggcttgctctctcttcacctgctctgttctacagcacactaccagaagacagcagaaatgaaaagcatttactttgtggctggattatttgtaatgctggtacaaggcagctggcaacgttcccttcaagacacagaggagaaatccagatcattctcagcttcccaggcagacccactcagtgatcctgatcagatgaacgaggacaagcgccattcacagggcacattcaccagtgactacagcaagtatctggactccaggcgtgcccaagattttgtgcagtggttgatgaataccaagaggaacaggaataacattgccaaacgtcacgatgaatttgagagacatgctgaagggacctttaccagtgatgtaagttcttatttggaaggccaagctgccaaggaattcattgcttggctggtgaaaggccgaggaaggcgagatttcccagaagaggtcgccattgttgaagaacttggccgcagacatgctgatggttctttctctgatgagatgaacaccattcttgataatcttgccgccagggactttataaactggttgattcagaccaaaatcactgacaggaaataactatatcactattcaagatcatcttcacaacatcacctgctagccacgtgggatgtttgaaatgttaagtcctgtaaatttaagaggtgtattctgaggccacattgctttgcatgccaataaataaattttcttttagtgttgtgtagccaaaaattacaaatggaataaagttttatcaaaatattgctaaaatatcagctttaaaatatgaaagtgctagattctgttattttcttcttattttggatgaagtaccccaacctgtttacatttagcgataaaattatttttctatgatataatttgtaaatgtaaattattccgatctgacatatctgcattataataataggagaatagaagaactggtagccacagtggtgaaattggaaagagaactttcttcctgaaacctttgtcttaaaaatactcagctttcaatgtatcaaagatacaattaaataaaattttcaagcttctttaccattgtct"
        ##gcg_t1 = Transcript(gcg_ts, "NM_002054.4", {344: gcg_v1})
        #gcg_t1 = Transcript(self.gcg_ts, 'GLUC_HUMAN', "NM_002054.4", [self.gcg_v1])
        #self.w_t = Peptide("TRANSCRIPT", {gcg_p1: [0]})

    def test_consistency(self):
        """
        tests all __*__ (including init)
        test has several asserts! If one fails, the following will not be evaluated!
        """
        self.assertTrue(repr(self.simple) == "PEPTIDE:\n SYFPEITHI")
        self.assertTrue(repr(self.w_p) == "PEPTIDE:\n PROTEIN\nin TRANSCRIPT: GLUC_HUMAN\n\tVARIANTS:\nin PROTEIN: GLUC_HUMAN")
        self.assertTrue(repr(self.w_v) == "PEPTIDE:\n VARIANT\nin TRANSCRIPT: GLUC_HUMAN\n\tVARIANTS:\n\tVariant(g.162145588G>T)\nin PROTEIN: GLUC_HUMAN")

    def test_getitem(self):
        self.assertTrue(self.simple[1:3] == 'YF')
        #TODO: document to have variant peptides from Protein with Variants use Generator

    def test_get_all_variants(self):
        self.assertTrue(repr(self.w_v.get_variants_by_protein("GLUC_HUMAN")) == repr([self.gcg_v1]))

    def test_get_all_proteins(self):
        self.assertTrue(repr(self.simple.get_all_proteins()) == repr([]))
        self.assertTrue(repr(self.w_p.get_all_proteins()) == repr([self.gcg_p1]))

    def test_get_all_transcripts(self):
        self.assertTrue(repr(self.w_v.get_all_transcripts()) == repr([Transcript(_seq="", _transcript_id="GLUC_HUMAN")]))
        self.assertTrue(repr(self.w_p.get_all_transcripts()) == repr([self.gcg_t1]))
