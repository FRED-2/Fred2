from unittest import TestCase
from Fred2.Core import Allele

__author__ = 'walzer'


class TestAllele(TestCase):
    def setUp(self):
        self.simple = Allele("HLA-A*02:01")

    def test_consistency(self):
        """
        tests all __*__ (including init)
        test has several asserts! If one fails, the following will not be evaluated!
        """
        self.assertTrue(repr(self.simple) == "HLA-A*02:01")
        self.assertEqual(self.simple, Allele("HLA-A*02:01"))
        self.assertNotEqual(repr(self.simple), Allele("HLA-A*02:01:666"))

