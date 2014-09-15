# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'schubert'

'''
Unit Test for Transcript Generator
'''
import unittest

from Fred2.Core.Variant import Variant, VariationType, MutationSyntax
from Fred2.Core.Generator import generate_transcripts_from_variants
from Fred2.IO.MartsAdapter import MartsAdapter
from Fred2.IO.FileReader import read_annovar_exonic

from Fred2.TestSuite.DummyAdapter import DummyAdapter
from Fred2.TestSuite.VariantsForTesting import *


class TranskriptGeneratorTestCase(unittest.TestCase):

    def setUp(self):
        self.trid = "NM_001114377" #FOXP3
        #self, id, type, chrom, genomePos, ref, obs, coding, isHomozygous, 
        # isSynonymous, metadata=None)
        self.non_syn_hetero_snp = Variant("COSM1122493",VariationType.SNP,"X",
                                          49111949, "G", "T",
                                          {"NM_001114377":MutationSyntax( \
                                            "NM_001114377", 757, 218, "", "")
                                          }, False, False)

        self.non_frame_shift_del = Variant("COSM1122495",VariationType.DEL,"X", 
                                           49113232, "CTT", "",
                                           {"NM_001114377":MutationSyntax( \
                                            "NM_001114377", 616, 206, "", "")
                                           }, True, False)

        self.syn_homo_snp = Variant("COSM1122494",VariationType.SNP,"X",
                                    49112257, "G", "A",
                                    {"NM_001114377":MutationSyntax( \
                                     "NM_001114377", 654, 218, "", "")
                                    }, False, True)

        self.db_adapter = MartsAdapter()


    # def test_non_syn_hetero_snp_trans_number(self):
    #     """
    #     tests if the number of generated transcripts for a heterozygous
    #     transcript is correct
    #
    #     1 hetero vars = 2 transcripts
    #     :return:
    #     """
    #     vars_ = \
    #     [self.non_syn_hetero_snp, self.non_frame_shift_del,self.syn_homo_snp]
    #
    #     trans = \
    #     [t for t in generate_transcripts_from_variants(vars_, self.db_adapter)]
    #     # print trans
    #
    #     self.assertTrue(len(trans) == 2**sum(not v.isHomozygous for v in vars_))
    #
    # def test_simple_incorporation(self):
    #     """
    #     test simple variant incorporation. only 1 variant in 1 transcript.
    #     input reference transcript: AAAAACCCCCGGGGG
    #
    #     variant 3: insert TT after pos 7
    #
    #     variant 1: SNP C -> T at pos 2
    #
    #     variant 4: del CCCCC after pos 9
    #     """
    #     dummy_db = DummyAdapter()
    #
    #     # INSERTIONS:
    #     dummy_vars = [ var_3]
    #     trans = generate_transcripts_from_variants(dummy_vars, dummy_db).next()
    #
    #     self.assertEqual(str(trans), "AAAAACCTTCCCGGGGG")
    #
    #     # SNPs:
    #     dummy_vars = [ var_1]
    #     trans = generate_transcripts_from_variants(dummy_vars, dummy_db).next()
    #     self.assertEqual(str(trans), "ATAAACCCCCGGGGG")
    #
    #
    #     # DELETIONS:
    #     dummy_vars = [ var_4]
    #     trans = generate_transcripts_from_variants(dummy_vars, dummy_db).next()
    #     self.assertEqual(str(trans), "AAAAACCCCG")
    #
    #
    #
    # def test_offset_single(self):
    #     """
    #     tests if offset is correctly handled when several variants for one
    #     transcript occur. still only one transcript with one transcript variant.
    #     reference transcript: AAAAACCCCCGGGGG
    #
    #     Each variant so that it is clearly down stream of
    #     it's predecessor
    #
    #     """
    #     dummy_db = DummyAdapter()
    #
    #     # 1) INS, SNP, DEL
    #     dummy_vars = [ var_3, var_7, var_6]
    #     trans = generate_transcripts_from_variants(dummy_vars, dummy_db).next()
    #
    #     self.assertEqual(str(trans), "AAAAACCTTCTCGGG")
    #
    #     # 2.) INS, DEL, INS
    #     dummy_vars = [ var_9, var_4, var_8]
    #     trans = generate_transcripts_from_variants(dummy_vars, dummy_db).next()
    #     self.assertEqual(str(trans), "AATTAAACCCCGTTT")
    #
    #
    # def test_heterozygous_variants(self):
    #     """
    #     Create multiple transcript variants for a transcript, given a set
    #     containing heterozygous variants .
    #
    #     Variants:
    #     3-DEL(-2)  , 5-INS(+3)  , 7-DEL(-4)
    #     HET-DEL(-2), HOM-INS(+3), HET-DEL(-1)
    #
    #     Reference sequence:
    #     AAAAACCCCCGGGGG
    #     AAATTTCGGGGG (DEL,INS,DEL)
    #     AAATTTCCCCCGGGGG (DEL,INS)
    #     AAAAATTTCCGGGG (INS,DEL)
    #     AAAAATTTCCCCCGGGGG (INS)
    #
    #     GGGGGCCCCCAAAAA
    #     GGGTTTCAAAAA (DEL,INS,DEL)
    #     GGGTTTCCCCCAAAAA (DEL,INS)
    #     GGGGGTTTCAAAAA (INS,DEL)
    #     GGGGGTTTCCCCCAAAAA (INS)
    #     """
    #
    #     dummy_db = DummyAdapter()
    #
    #     # 1) INS, SNP, DEL
    #     dummy_vars = [var_10, var_11, var_12]
    #     trans_gener = generate_transcripts_from_variants(dummy_vars, dummy_db)
    #     trans = [t for t in trans_gener]
    #
    #
    #     trans = map(str, trans)
    #
    #     self.assertEqual(len(trans), 8)
    #     self.assertTrue("AAATTTCCGGGG" in trans)
    #     self.assertTrue("AAATTTCCCCCGGGGG" in trans)
    #     self.assertTrue("AAAAATTTCCGGGG" in trans)
    #     self.assertTrue("AAAAATTTCCCCCGGGGG" in trans)
    #
    #     self.assertTrue("GGGTTTCCAAAA" in trans)
    #     self.assertTrue("GGGTTTCCAAAA" in trans)
    #     self.assertTrue("GGGTTTCCCCCAAAAA" in trans)
    #     self.assertTrue("GGGGGTTTCCCCCAAAAA" in trans)

    def test_varinat_reader(self):
        vars = read_annovar_exonic("../IO/snp_annot_donor.txt.exonic_variant_function", gene_filter=[])
        print vars, len(vars)
        trans = list(generate_transcripts_from_variants(vars, self.db_adapter))
        print trans
        transToVar = {}
        for v in vars:
            for trans_id in v.coding.iterkeys():
                transToVar.setdefault(trans_id, []).append(v)


if __name__ == '__main__':
    unittest.main()

