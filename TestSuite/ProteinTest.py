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

    def test1_protein_construction_novariants(self):
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


    def test2_generate_peptides_novariants(self):
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
    def test3_protein_from_variants(self):
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


    def test4_peptides_from_variants(self):
        """
        Ref trancript: AAAAACCCCCGGGGG
        ref protein:   KNPRG
        ref peps(3):   KNPR, NPRG

        variant1: heterozygous, fs+1 in first aa
        variant2: heterozygous, insertion +2 in last aa

        trans-var1: TKPPGA
        1: peps(3): TKPP, KPPG, PPGA

        trans-var2: KNPRG
        2: peps(3): KNPR, NPRG

        Output:
        -------
        PEPTIDE: PPGA
            TRANSCRIPT: tsc_1:FRED2_3
                 Variant(15CC)
                 Variant(1C)
        PEPTIDE: KPPG
            TRANSCRIPT: tsc_1:FRED2_3
                 Variant(1C)
        PEPTIDE: TKPP
            TRANSCRIPT: tsc_1:FRED2_3
                 Variant(1C)
        
        PEPTIDE: KNPR
            TRANSCRIPT: tsc_1:FRED2_0
        PEPTIDE: NPRG
            TRANSCRIPT: tsc_1:FRED2_0
        """

        peps_trans1 = ["KNPR", "NPRG"]
        peps_trans2 = ["PPGA", "KPPG", "TKPP"]
        expected_vars = ["Variant(1C)", "Variant(15CC)"]
        expected = peps_trans1 + peps_trans2

        dummy_db = DummyAdapter()
        dummy_vars = [var_13, var_14]

        proteins = []
        for trans in generate_transcripts_from_variants(dummy_vars, dummy_db):
            ### GET PROTS:
            # IGNORE invalid sequence lengths
            try:
                proteins.append(trans.translate())
            except ValueError:
                pass

        peptides = generate_peptides_from_protein(proteins, 4)

        sequences = [str(pep) for pep in peptides]

        # Check if all peptides are generated as expected
        self.assertTrue(all(pep in sequences for pep in expected))
        # no duplicates or more than the expected ones:
        self.assertEqual(len(peptides), len(expected))

        vari_peps = [pep.get_all_variants() for pep in peptides \
                     if str(pep) in peps_trans2]

        vars_ = [str(var) for varlist in vari_peps for var in varlist]

        # Check that for the peptides from the transcript containing the 
        # variants, we also get all expected variants. Especally the first
        # variant needs to be present in all peptides
        self.assertTrue(all(var in vars_ for var in expected_vars))


if __name__ == '__main__':
    unittest.main()
