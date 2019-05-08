"""Unit test for protein class
"""

import unittest

# Variants and Generator
from Fred2.Core.Generator import generate_transcripts_from_variants
from Fred2.test.DummyAdapter import DummyAdapter
from Fred2.test.VariantsForTesting import *

from Fred2.Core.Protein import Protein
from Fred2.Core.Generator import generate_peptides_from_proteins
from Fred2.Core.Generator import generate_proteins_from_transcripts
from Fred2.IO.ADBAdapter import EIdentifierTypes

class TestProteinClass(unittest.TestCase):
    def setUp(self):
        # generate a Protein for use in the tests of this class
        self.single_protein = Protein("ASDERWQTGHKILPMNVFCY", gene_id='gene 1', transcript_id='someID')

        # generate a set of Proteins for use in the tests of this class
        self.prot_set = list()
        self.prot_set.append(Protein("IIIVRC", gene_id='gene 1', transcript_id='set entry 1'))
        self.prot_set.append(Protein("VRCVR", gene_id='gene 1', transcript_id='set entry 2'))
        self.prot_set.append(Protein("IIVRCIT", gene_id='gene 1', transcript_id='set entry 3'))
        self.prot_set.append(Protein("IVRC", gene_id='gene 1', transcript_id='set entry 4'))

    def test1_protein_construction_novariants(self):
        """
        Test if a single protein is correctly constructed from a given sequence and gene ID.
        """
        self.assertEqual(self.single_protein.gene_id, "gene 1", 'incorrect Id')
        self.assertEqual(self.single_protein.transcript_id, "someID", 'incorrect transcript id')
        self.assertEqual(str(self.single_protein), "ASDERWQTGHKILPMNVFCY", 'incorrect sequence')
        self.assertIsNone(self.single_protein.orig_transcript)
        self.assertEqual(self.single_protein.vars, {})

    def test2_generate_peptides_novariants(self):
        """
        Test if a list of proteins is correctly broken into peptide fragments.
        Here the proteins are constructed just from their sequence, having no
        transcript or variant information.
        """
        pep_set = generate_peptides_from_proteins(self.prot_set, 3)

        # # Print peptide generator results:
        # for pep in pep_set:
        #     print pep, pep.proteins.items()
        #     print pep, pep.vars.items()
        #     print pep, pep.transcripts.items()

        # get the number of peptides generated for each protein in self.prot_set and sum up
        number_of_peps = sum(len(list(pep.proteins.keys())) for pep in pep_set)
        # The total number of peptides of length 3 from all proteins in self.pro_set should be 14
        self.assertEqual(number_of_peps, 14)

        # generated pep_set should consist only of unique-sequence entries
        unique_test_prot_set = list()
        unique_test_prot_set.extend(self.prot_set)
        unique_test_prot_set.extend(self.prot_set)

        unique_test_pep_set = set(generate_peptides_from_proteins(unique_test_prot_set, 3))
        unique_test_pep_seqs = set([str(pep) for pep in unique_test_pep_set])
        self.assertEqual(len(unique_test_pep_set), len(unique_test_pep_seqs))

    def test3_protein_from_variants(self):
        """
        Generate some transcripts from the 3 input variants
        (should give 8 transcripts, check also if all fields are complete)
        Using a protein made from variants:

        Translate to proteins (check if all fields are there/filled)

        fragment to unique peptides
        (check for uniqueness of sequences, check fields of peptides, check
        correctness of fragments)
        """
        dummy_db = DummyAdapter()
        dummy_vars = [var_10, var_11, var_12]

        proteins = []
        t = list(generate_transcripts_from_variants(dummy_vars, dummy_db, EIdentifierTypes.REFSEQ))
        for trans in t:
            # check gene id field:
            print(trans)
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
                proteins.append(next(generate_proteins_from_transcripts(trans)))
            except ValueError:
                pass

        self.assertEqual(len(proteins), 8)

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
            self.assertEqual(len(set(e for subl in prot.vars.values() for e in subl)), len(orig.vars))

            # check sequence:
            self.assertTrue(str(prot) > 2)

        ## GENERATE Peptides:
        peptides = generate_peptides_from_proteins(proteins,2)


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
        #TODO Somewhere here a print statement is called
        peps_trans1 = ["KNPR", "NPRG"]
        peps_trans2 = ["PPGA", "KPPG", "TKPP"]
        expected_vars = ["Variant(1C)", "Variant(15CC)"]
        expected = peps_trans1 + peps_trans2

        dummy_db = DummyAdapter()
        dummy_vars = [var_13, var_14]

        proteins = []
        transcripts = list(generate_transcripts_from_variants(dummy_vars, dummy_db, EIdentifierTypes.REFSEQ))
        for trans in transcripts:
            ### GET PROTS:
            # IGNORE invalid sequence lengths
            try:
                proteins.append(next(generate_proteins_from_transcripts(trans)))
            except ValueError:
                pass

        peptides = list(generate_peptides_from_proteins(proteins, 4))

        sequences = [str(pep) for pep in peptides]

        # Check if all peptides are generated as expected
        self.assertTrue(all(pep in sequences for pep in expected))
        # no duplicates or more than the expected ones:
        self.assertEqual(len(peptides), len(expected))

        #vari_peps = [pep.get_all_variants() for pep in peptides \
        #             if str(pep) in peps_trans2]

        #vars_ = [str(var) for varlist in vari_peps for var in varlist]

        # Check that for the peptides from the transcript containing the
        # variants, we also get all expected variants. Especally the first
        # variant needs to be present in all peptides
        for prot in proteins:
            for p in peptides:
                try:
                    vars_ = list(map(str, p.get_variants_by_protein(prot.transcript_id)))
                    expected_vars = [str(v) for vars in prot.vars.values() for v in vars]
                    print("peptide vars: ", vars_)
                    print("Prot vars: ", expected_vars)
                    print(repr(p))
                    print(repr(prot))
                    self.assertTrue(all(var in expected_vars for var in vars_))
                except KeyError:
                    pass


if __name__ == '__main__':
    unittest.main()
