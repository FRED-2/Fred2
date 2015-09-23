from unittest import TestCase

__author__ = 'walzer'

from Fred2.Core import Variant
from Fred2.Core import VariationType
from Fred2.Core import MutationSyntax


class TestVariant(TestCase):
    def setUp(self):
        self.simple = Variant("rs5650", VariationType.SNP, 2, 162145588, 'G', 'T',
                        {"NM_002054.4": MutationSyntax("NM_002054.4", 344, 115, "c.344C>A", "p.A115D")}, False, False)

    def test_consistency(self):
        """
        tests all __*__ (including init)
        test has several asserts! If one fails, the following will not be evaluated!
        """
        self.assertTrue(repr(self.simple) == "Variant(g.162145588G>T)")

    def test_get_transcript_offset(self):
        self.assertTrue(self.simple.get_transcript_offset() == 0)

    def test_get_shift(self):
        self.assertTrue(self.simple.get_transcript_offset() == 0)

    def test_get_transcript_position(self):
        self.assertTrue(self.simple.get_annotated_transcript_pos("NM_002054.4") == 344)

    def test_get_protein_position(self):
        self.assertTrue(self.simple.get_annotated_protein_pos("NM_002054.4") == 115)