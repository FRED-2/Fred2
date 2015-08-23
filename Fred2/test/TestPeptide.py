from unittest import TestCase

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
        gcg_p1 = Protein(self.gcg_ps, 'GLUC_HUMAN')
        self.w_p = Peptide("PROTEIN", {'GCG': gcg_p1})

        gcg_v1 = Variant("rs5650", VariationType.SNP, 2, 162145588, 'G', 'T',
                         {"NM_002054.4": MutationSyntax("NM_002054.4", 344, 115, "c.344C>A", "p.A115D")}, False, False)
        self.w_v = Peptide("VARIANT", {'GCG': gcg_p1}, {'NM_002054.4': [gcg_v1]})

        self.gcg_ts = "gcatagaatgcagatgagcaaagtgagtgggagagggaagtcatttgtaacaaaaactcattatttacagatgagaaatttatattgtcagcgtaatatctgtgaggctaaacagagctggagagtatataaaagcagtgcgccttggtgcagaagtacagagcttaggacacagagcacatcaaaagttcccaaagagggcttgctctctcttcacctgctctgttctacagcacactaccagaagacagcagaaatgaaaagcatttactttgtggctggattatttgtaatgctggtacaaggcagctggcaacgttcccttcaagacacagaggagaaatccagatcattctcagcttcccaggcagacccactcagtgatcctgatcagatgaacgaggacaagcgccattcacagggcacattcaccagtgactacagcaagtatctggactccaggcgtgcccaagattttgtgcagtggttgatgaataccaagaggaacaggaataacattgccaaacgtcacgatgaatttgagagacatgctgaagggacctttaccagtgatgtaagttcttatttggaaggccaagctgccaaggaattcattgcttggctggtgaaaggccgaggaaggcgagatttcccagaagaggtcgccattgttgaagaacttggccgcagacatgctgatggttctttctctgatgagatgaacaccattcttgataatcttgccgccagggactttataaactggttgattcagaccaaaatcactgacaggaaataactatatcactattcaagatcatcttcacaacatcacctgctagccacgtgggatgtttgaaatgttaagtcctgtaaatttaagaggtgtattctgaggccacattgctttgcatgccaataaataaattttcttttagtgttgtgtagccaaaaattacaaatggaataaagttttatcaaaatattgctaaaatatcagctttaaaatatgaaagtgctagattctgttattttcttcttattttggatgaagtaccccaacctgtttacatttagcgataaaattatttttctatgatataatttgtaaatgtaaattattccgatctgacatatctgcattataataataggagaatagaagaactggtagccacagtggtgaaattggaaagagaactttcttcctgaaacctttgtcttaaaaatactcagctttcaatgtatcaaagatacaattaaataaaattttcaagcttctttaccattgtct"
        #gcg_t1 = Transcript(gcg_ts, "NM_002054.4", {344: gcg_v1})
        gcg_t1 = Transcript(self.gcg_ts, 'GLUC_HUMAN', "NM_002054.4", [gcg_v1])
        self.w_t = Peptide("TRANSCRIPT", {'GCG': gcg_p1}, {'NM_002054.4': [gcg_v1]}, {"NM_002054.4": gcg_t1})

    def test_consistency(self):
        """
        tests all __*__ (including init)
        test has several asserts! If one fails, the following will not be evaluated!
        """
        self.assertTrue(repr(self.simple) == "PEPTIDE:\n SYFPEITHI")
        self.assertTrue(repr(self.w_p) == "PEPTIDE:\n PROTEIN\nin PROTEIN: GCG")
        self.assertTrue(repr(self.w_v) == "PEPTIDE:\n VARIANT\nin TRANSCRIPT: NM_002054.4\n\tVARIANTS:\n\tVariant(g.162145588G>T)\nin PROTEIN: GCG")
        self.assertTrue(repr(self.w_t) == "PEPTIDE:\n TRANSCRIPT\nin TRANSCRIPT: NM_002054.4\n\tVARIANTS:\n\tVariant(g.162145588G>T)\nin PROTEIN: GCG")

    def test_getitem(self):
        self.assertTrue(self.simple[1:3] == 'YF')
        #TODO: document to have variant peptides from Protein with Variants use Generator

    def test_get_all_variants(self):
        gcg_v_test = Variant("rs5650", VariationType.SNP, 2, 162145588, 'G', 'T',
                         {"NM_002054.4": MutationSyntax("NM_002054.4", 344, 115, "c.344C>A", "p.A115D")}, False, False)
        self.assertTrue(repr(self.w_v.get_all_variants()) == repr([gcg_v_test]))

    def test_get_all_proteins(self):
        gcg_p_test = Protein(self.gcg_ps, 'GLUC_HUMAN', _transcript_id="Protein_1")

        self.assertTrue(repr(self.simple.get_all_proteins()) == repr([]))
        self.assertTrue(repr(self.w_p.get_all_proteins()) == repr([gcg_p_test]))

    def test_get_all_transcripts(self):
        gcg_v_test = Variant("rs5650", VariationType.SNP, 2, 162145588, 'G', 'T',
                         {"NM_002054.4": MutationSyntax("NM_002054.4", 344, 115, "c.344C>A", "p.A115D")}, False, False)
        gcg_t_test = Transcript(self.gcg_ts, 'GLUC_HUMAN', "NM_002054.4", [gcg_v_test])

        self.assertTrue(repr(self.w_v.get_all_transcripts()) == repr([Transcript(_seq="", _transcript_id="NM_002054.4")]))
        self.assertTrue(repr(self.w_t.get_all_transcripts()) == repr([gcg_t_test]))
