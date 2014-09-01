"""Unit test for protein class
"""

import unittest

# Variants and Generator
from Fred2.Core.Generator import generate_transcripts_from_variants
from Fred2.IO.DummyAdapter import DummyAdapter
from Fred2.UnitTest.VariantsForTesting import *

# Transcripts

# Protein
from Fred2.Core.Protein import Protein
from Fred2.Core.Generator import generate_peptides_from_protein

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
        self.assertEqual(self.single_protein.vars, {})


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
        #     print pep, pep.proteins.items()
        #     print pep, pep.vars.items()
        #     print pep, pep.transcripts.items()

        # The total number of fragments should be 14
        # which is the sum over the individual originating proteins
        self.assertEqual(get_total_peps(pep_set), 14)

        # pep_set consists only of unique entries
        pep_unique_seq = set([str(pep) for pep in pep_set])
        self.assertEqual(len(pep_set), len(pep_unique_seq))



    # Using a protein made from variants:
    def test_2_protein_from_variants(self):
        """
        Generate some transcripts from the 3 input variants
        (should give 8 transcripts, check also if all fields are complete)

        Translate to proteins (check if all fields are there/filled)

        fragment to unique peptides
        (check for uniqueness of sequences, check fields of peptides, check 
        correctness of fragments)
        """
        dummy_db = DummyAdapter()
        dummy_vars = [var_10, var_11, var_12]

        proteins = []
        for trans in generate_transcripts_from_variants(dummy_vars, dummy_db):
            # check gene id field:
            self.assertEqual(trans.gene_id, "gene_1")

            # check trans id name:
            name = trans.transcript_id.split(":FRED2_")
            self.assertEqual(len(name), 2)
            self.assertTrue(name[0] == "tsc_1" or name[0] == "tsc_2")
            self.assertTrue(len(name[1]) == 1 and name[1].isdigit)

            # check var:
            self.assertIsNotNone(trans.vars)
            self.assertTrue(len(trans.vars) > 0)

            # check sequence:
            self.assertTrue(str(trans) > 5)

            ### GET PROTS:
            # IGNORE invalid sequence lengths
            try:
                proteins.append(trans.translate())
            except ValueError:
                pass
                #print "invalid transcripts, due to non multiple of 3 length: ",\
                #      trans

        self.assertEqual(len(proteins), 4)

        ## CHECK Proteins:
        for prot in proteins:
            self.assertEqual(prot.gene_id, "gene_1")

            # check trans id name:
            name = prot.transcript_id.split(":FRED2_")
            self.assertEqual(len(name), 2)
            self.assertTrue(name[0] == "tsc_1" or name[0] == "tsc_2")
            self.assertTrue(len(name[1]) == 1 and name[1].isdigit)

            orig = prot.orig_transcript
            self.assertEqual(prot.transcript_id, orig.transcript_id)

            self.assertEqual(len(prot.vars), len(orig.vars))

            # check sequence:
            self.assertTrue(str(prot) > 2)

        ## GENERATE Peptides:
        peptides = generate_peptides_from_protein(proteins,2)

        print peptides
        print len(peptides)




if __name__ == '__main__':
    unittest.main()
