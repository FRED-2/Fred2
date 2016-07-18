from unittest import TestCase
from Fred2.Core import Allele, CombinedAllele

__author__ = 'walzer,schubert'


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

    def test_combined_allele(self):
        comb = CombinedAllele("HLA-DPA1*01:03-DPB1*01:01")
        self.assertTrue(repr(comb) == "HLA-DPA1*01:03-DPB1*01:01")
        self.assertEqual(comb, CombinedAllele("HLA-DPA1*01:03-DPB1*01:01"))
        self.assertNotEqual(comb, CombinedAllele("HLA-DPA1*01:03-DPB1*01:02"))

    def test_allele_factory(self):
        a = Allele("HLA-DPA1*01:03-DPB1*01:01", prob=1)
        b = Allele("HLA-A*02:01", prob=2)
        self.assertIsInstance(a, CombinedAllele)
        self.assertEqual(a.prob,1)
        self.assertIsInstance(b, Allele)
        self.assertEqual(b.prob, 2)