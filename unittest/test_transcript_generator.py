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


class TranskriptGeneratorTestCase(unittest.TestCase):

    def setUp(self):
        self.trid = "NM_001114377" #FOXP3
        #self, id, type, chrom, genomePos, ref, obs, coding, isHomozygous, isSynonymous, metadata=None)
        self.non_syn_hetero_snp = Variant("COSM1122493",VariationType.SNP,"X", 49111949, "G", "T",
                                     {"NM_001114377":MutationSyntax("NM_001114377",840,218,"","")}, False, False)

        self.non_frame_shift_del = Variant("COSM1122495",VariationType.DEL,"X", 49113237, "CTT", "",
                                      {"NM_001114377":MutationSyntax("NM_001114377",699,171,"","")}, True, False)

        self.syn_homo_snp = Variant("COSM1122494",VariationType.SNP,"X", 49112257, "G", "A",
                               {"NM_001114377":MutationSyntax("NM_001114377",737, 183,"","")}, True, True)

        self.db_adapter = MartsAdapter()

    def test_non_syn_hetero_snp_trans_number(self):
        """
        tests if the number of generated transcripts for a heterozygous transcript is correct

        1 hetero vars = 2 transcripts
        :return:
        """
        vars = [self.non_syn_hetero_snp]
        trans = [t for t in generate_transcripts_from_variants(vars, self.db_adapter)]
        print trans
        self.assertTrue(len(trans) == 2)



if __name__ == '__main__':
    unittest.main()

