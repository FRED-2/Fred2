# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'schubert'

'''
Unit Test for Transcript Generator
'''
import unittest

from Fred2.Core.Base import COMPLEMENT
from Fred2.Core.Variant import Variant, VariationType, MutationSyntax
from Fred2.Core.Generator import generate_transcripts_from_variants,_incorp
from Fred2.IO.MartsAdapter import MartsAdapter
from Fred2.IO.ADBAdapter import EAdapterFields

from Fred2.IO.DummyAdapter import DummyAdapter
from Fred2.TestSuite.VariantsForTesting import *

class TranskriptGeneratorTestCase(unittest.TestCase):

    def setUp(self):
        self.trid = "NM_001114377" #FOXP3
        #self, id, type, chrom, genomePos, ref, obs, coding, isHomozygous, isSynonymous, metadata=None)
        self.non_syn_hetero_snp = Variant("COSM1122493",VariationType.SNP,"X", 49111949, "G", "T",
                                     {"NM_001114377":MutationSyntax("NM_001114377",757,218,"","")}, False, False)

        self.non_frame_shift_del = Variant("COSM1122495",VariationType.DEL,"X", 49113232, "CTT", "",
                                      {"NM_001114377":MutationSyntax("NM_001114377",616,206,"","")}, True, False)

        self.syn_homo_snp = Variant("COSM1122494",VariationType.SNP,"X", 49112257, "G", "A",
                               {"NM_001114377":MutationSyntax("NM_001114377",654, 218,"","")}, False, True)

        self.db_adapter = MartsAdapter()


    def test_non_syn_hetero_snp_trans_number(self):
        """
        tests if the number of generated transcripts for a heterozygous transcript is correct

        1 hetero vars = 2 transcripts
        :return:
        """
        vars = [self.non_syn_hetero_snp, self.non_frame_shift_del,self.syn_homo_snp]
        trans = [t for t in generate_transcripts_from_variants(vars, self.db_adapter)]
        # print trans
        # print "\ndiff: ", [ (a,b) for a,b in zip(str(trans[0]), str(trans[1])) if a != b]
        self.assertTrue(len(trans) == 2**sum(not v.isHomozygous for v in vars))


    # def test_trans_post_in_var_is_correct(self):
    #     """
    #     :return:
    #     """
    #     vars = [self.syn_homo_snp]
    #
    #     query = self.db_adapter.get_transcript_information(self.trid)
    #     trans = query[EAdapterFields.SEQ]
    #     geneid = query[EAdapterFields.GENE]
    #     strand = query[EAdapterFields.STRAND]
    #
    #     #if its a reverse transcript form the complement of the variants
    #     if strand == "-":
    #         for v in vars:
    #             v.ref = v.ref[::-1].translate(COMPLEMENT)
    #             v.obs = v.obs[::-1].translate(COMPLEMENT)
    #     print
    #     for v in vars:
    #         print "SNP reference %s at transcript pos %i , transcript base %s at same position"% (v.ref,
    #                                         v.get_transcript_position(self.trid),
    #                                         trans[v.get_transcript_position(self.trid)])
    #         print
    #         self.assertTrue(v.ref == trans[v.get_transcript_position(self.trid)-1])

    def test_simple_incorporation(self):
        """
        test simple variant incorporation. only 1 variant in 1 transcript.
        input reference transcript: AAAAACCCCCGGGGG

        variant 3: insert TT after pos 7

        variant 1: SNP C -> T at pos 2

        variant 4: del CCCCC after pos 9
        """
        dummy_db = DummyAdapter()

        # INSERTIONS:
        dummy_vars = [ var_3]
        trans = generate_transcripts_from_variants(dummy_vars, dummy_db).next()
        
        self.assertEqual(str(trans), "AAAAACCTTCCCGGGGG")

        # SNPs:
        dummy_vars = [ var_1]
        trans = generate_transcripts_from_variants(dummy_vars, dummy_db).next()
        self.assertEqual(str(trans), "ATAAACCCCCGGGGG")


        # DELETIONS:
        dummy_vars = [ var_4]
        trans = generate_transcripts_from_variants(dummy_vars, dummy_db).next()
        self.assertEqual(str(trans), "AAAAACCCCG")



    def test_offset_single(self):
        """
        tests if offset is correctly handled when several variants for one 
        transcript occur. still only one transcript with one transcript variant.
        reference transcript: AAAAACCCCCGGGGG

        Each variant so that it is clearly down stream of
        it's predecessor

        """
        dummy_db = DummyAdapter()

        # 1) INS, SNP, DEL
        dummy_vars = [ var_3, var_7, var_6]
        trans = generate_transcripts_from_variants(dummy_vars, dummy_db).next()
        
        self.assertEqual(str(trans), "AAAAACCTTCTCGGG")

        # 2.) INS, DEL, INS
        dummy_vars = [ var_9, var_4, var_8]
        trans = generate_transcripts_from_variants(dummy_vars, dummy_db).next()
        self.assertEqual(str(trans), "AATTAAACCCCGTTT")


    def test_heterozygous_variants(self):
        """
        Create multiple transcript variants for a transcript, given a set
        containing heterozygous variants .

        Variants:
        3-DEL(-2)  , 5-INS(+3)  , 7-DEL(-4)
        HET-DEL(-2), HOM-INS(+3), HET-DEL(-1)

        Reference sequence:
        AAAAACCCCCGGGGG
        AAATTTCGGGGG (DEL,INS,DEL)
        AAATTTCCCCCGGGGG (DEL,INS)
        AAAAATTTCGGGGG (INS,DEL)
        AAAAATTTCCCCCGGGGG (INS)

        GGGGGCCCCCAAAAA
        GGGTTTCAAAAA (DEL,INS,DEL)
        GGGTTTCCCCCAAAAA (DEL,INS)
        GGGGGTTTCAAAAA (INS,DEL)
        GGGGGTTTCCCCCAAAAA (INS)
        """

        dummy_db = DummyAdapter()

        # 1) INS, SNP, DEL
        dummy_vars = [var_10, var_11, var_12]
        trans_gener = generate_transcripts_from_variants(dummy_vars, dummy_db)
        trans = [t for t in trans_gener]

        print "\n\n\nTEST HETERO\n\n\n", trans

        trans = map(str, trans)

        self.assertEqual(len(trans), 8)
        self.assertTrue("AAATTTCGGGGG" in trans)
        self.assertTrue("AAATTTCCCCCGGGGG" in trans)
        self.assertTrue("AAAAATTTCGGGGG" in trans)
        self.assertTrue("AAAAATTTCCCCCGGGGG" in trans)

        self.assertTrue("GGGTTTCAAAAA" in trans)
        self.assertTrue("GGGTTTCAAAAA" in trans)
        self.assertTrue("GGGTTTCCCCCAAAAA" in trans)
        self.assertTrue("GGGGGTTTCCCCCAAAAA" in trans)

if __name__ == '__main__':
    unittest.main()

