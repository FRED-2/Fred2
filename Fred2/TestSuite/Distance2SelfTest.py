"""
    Unittest for Distance2Self computation
"""
__author__ = 'mohr'


import unittest
from Fred2.Distance2Self.Distance2Self import Distance2Self
from Fred2.Distance2Self.DistanceMatrix import DistanceMatrix
from Fred2.Core.Peptide import Peptide
from Fred2.Data import DistanceMatrices



class Distance2SelfTest(unittest.TestCase):

    def setUp(self):
        self.peptides = [Peptide("SYFPEITHI"),Peptide("IHTIEPFYS")]
        self.testFasta = open("../testSequences.fasta", 'r')

    def test_init_matrices(self):
        blosum45 = DistanceMatrix(DistanceMatrices.BLOSUM45_distances)
        blosum50 = DistanceMatrix(DistanceMatrices.BLOSUM50_distances)
        blosum90 = DistanceMatrix(DistanceMatrices.BLOSUM90_distances)

        print blosum45
        print blosum50
        print blosum90

    def test_init(self):
        d2s = Distance2Self()

    def test_trie_generation(self):
        d2s = Distance2Self()
        d2s.generate_trie(self.testFasta)

    def test_distance_computation(self):
        d2s=Distance2Self()
        d2s.generate_trie(self.testFasta)
        result = d2s.calculate_distances(self.peptides, n=20)
        print result

        self.assertTrue(len(result.columns) == 22)
        self.assertTrue(len(result.index) == 2)


if __name__ == '__main__':
    unittest.main()