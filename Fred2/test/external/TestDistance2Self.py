"""
    Unittest for Distance2Self computation
"""
from Fred2.Core.Generator import generate_peptides_from_proteins

__author__ = 'mohr,walzer'


import unittest
from Fred2.Distance2Self.Distance2Self import Distance2Self
from Fred2.Distance2Self.DistanceMatrix import DistanceMatrix
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Protein import Protein
from Fred2.Data import DistanceMatrices
from Bio import SeqIO
import pkg_resources
from os import path


class Distance2SelfTest(unittest.TestCase):

    def setUp(self):
        self.peptides = [Peptide("SYFPEITHI"), Peptide("IHTIEPFYS")]
        testsequences_file = pkg_resources.resource_filename('Fred2', path.join('Data', 'examples', 'testSequences.fasta'))
        with open(testsequences_file, "rU") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        prot_set = [Protein(str(r.seq)) for r in records]
        unique_test_pep_set = generate_peptides_from_proteins(prot_set, 9)
        self.selfpeptides = [str(x) for x in unique_test_pep_set]

        small_prot_set = [Protein("MKERRIDMKEKNVKAKAPNKKVLGLTTKIFIALLAGAILGIVLCYLVPDSSFKKDVIVEGILYVIGQGFIRLMKMLVVPLVFCSLVCGSMAIGDTKKLGTVGVRTLAFYLATTALAVVALGVGNLINPGVGLDMSAIQSSAASVETMEATSLTDTILNIIPDNPINSLASGSMLQVIVFALIVGVILAKMGERAETVANFFSQFNDIMMEMTMMIMSLAPIGVFCLISRTFANIGFSAFIPLAKYMIGVLLALAIQCFGVYQILLKIFTGLNPIRFIKKFFPVMAFAFSTATSNATIPMSIDTLSKKVGVSKKISSFTIPLGATINMDGTSIMQGVAVVFAAQAFGIHLTPMDYVTVIGTATLASVGTAGVPSVGLVTLTMVFNSVGLPVEAIGLIMGIDRILDMTRTAVNITGDAVCTTIVAHQNGALDKKVFNETE"), Protein("MLKVWIAGASGQIGRALNDVLDPMQIEALNTDLDELDITDTDEVINFGTVNRPDVIINCTGITDTDECEANPEHAYRVNALGARNLSIVARKCGSKIVQLSTDDVFDGQSKKPYTEFDDTNPLTVYGRSKRAGENYVKEFTHKHFVIRSNWVYGHGGHNFVNRVLAAAEAGNGLSVASDQFGSPTSAKDLAKMIMYLISTNEYGTYHVTCRGVCSRYEFAQEILKLAGKDIELRAVPTEQSDLSAVRPPYAVLDNFILRIIEVYDMPDWKESLKEYMDERTED")]
        small_unique_test_pep_set = generate_peptides_from_proteins(small_prot_set, 9)
        self.fewselfpeptides = [str(x) for x in small_unique_test_pep_set]

    def test_init_matrices(self):
        blosum45 = DistanceMatrix(DistanceMatrices.BLOSUM45_distances, "BLOSUM45")
        blosum50 = DistanceMatrix(DistanceMatrices.BLOSUM50_distances, "BLOSUM50")
        blosum90 = DistanceMatrix(DistanceMatrices.BLOSUM90_distances, "BLOSUM90")

        print(blosum45)
        print(blosum50)
        print(blosum90)

    def test_init(self):
        blosum90 = DistanceMatrix(DistanceMatrices.BLOSUM90_distances, "BLOSUM90")
        d2s = Distance2Self(blosum90, self.selfpeptides)

    def test_distance_computation(self):
        import sys
        print("###", sys.getrecursionlimit())
        sys.setrecursionlimit(500000000)

        blosum90 = DistanceMatrix(DistanceMatrices.BLOSUM90_distances, "BLOSUM90")
        d2s=Distance2Self(blosum90, self.fewselfpeptides, 11)
        result = d2s.calculate_distances(self.peptides)
        print(result)

        self.assertTrue(len(result.columns) == 3)
        self.assertTrue(len(result.index) == 11)


if __name__ == '__main__':
    unittest.main()