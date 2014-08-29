"""Unit test for protein class
"""

from collections import OrderedDict

import unittest

from Fred2.Core.Protein import Protein, generate_peptides_from_protein

class TestProteinClass(unittest.TestCase):
    def setUp(self):
        # generate a Protein to test it
        self.single_protein = Protein("ASDERWQTGHKILPMNVFCY", 'gene 1', 
                                      'someID')

        self.prot_set = list()
        self.prot_set.append(Protein("IIIVRC", 'gene 1', 'set entry 1'))
        self.prot_set.append(Protein("VRCVR", 'gene 1', 'set entry 2'))
        self.prot_set.append(Protein("IIVRCIT", 'gene 1', 'set entry 3'))
        self.prot_set.append(Protein("IVRC", 'gene 1', 'set entry 4'))

    def test1_1_protein_construction_novariants(self):
        """
        Test if a single protein is correctly constructed from a given sequence
        and gene ID.
        (single_protein generation in 'setUp' method)
        """
        self.assertEqual(self.single_protein.gene_id, "gene 1", 'incorrect Id')
        self.assertEqual(self.single_protein.transcript_id, "someID", 
                         'incorrect transcript id')
        self.assertEqual(str(self.single_protein), "ASDERWQTGHKILPMNVFCY",
                         'incorrect sequence')
        self.assertIsNone(self.single_protein.orig_transcript)
        self.assertEqual(self.single_protein.vars, OrderedDict())


    def test1_2_generate_peptides_novariants(self):
        """
        Test if a list of proteins is correctly broken into peptide fragments.
        Here the proteins are constructed from their sequence, having no
        transcript or variant information.
        (prot_set generation in 'setUp' method)
        """
        def get_total_peps(pep_set):
            _sum = 0
            for pep in pep_set:
                _sum += len(pep.proteins.keys())
            return _sum

        pep_set = generate_peptides_from_protein(self.prot_set, 3)

        # # Print peptide generator results:
        # for pep in pep_set:
        #     print pep, pep.proteins.keys()
        for pep in pep_set:
            print pep, pep.proteins.iteritems()
            print pep, pep.vars.iteritems()
            print pep, pep.transcripts.iteritems()

        # The total number of fragments should be 14
        # which is the sum over the individual originating proteins
        self.assertEqual(get_total_peps(pep_set), 14)

        # pep_set consists only of unique entries
        pep_seqs = [str(pep) for pep in pep_set]
        self.assertEqual(len(pep_set), len(pep_seqs))

    # Using a protein made from transcripts:
    def test_2_1_protein_construction(self):
        pass

    def test_2_2_generate_peptides(self):
        pass


if __name__ == '__main__':
    unittest.main()
