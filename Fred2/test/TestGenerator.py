from unittest import TestCase

from Fred2.IO import MartsAdapter, read_annovar_exonic
from Fred2.Core.Variant import VariationType
from Fred2.test.DummyAdapter import DummyAdapter
from Fred2.test.VariantsForTesting import *
from Fred2.Core import Generator
from Fred2.IO.ADBAdapter import EIdentifierTypes
import os
import inspect
import Fred2


class GeneratorTest(TestCase):
    def setUp(self):
        self.trid = "NM_001114377"  # FOXP3
        # self, id, type, chrom, genomePos, ref, obs, coding, isHomozygous,
        # isSynonymous, metadata=None)
        self.non_syn_hetero_snp = Variant("COSM1122493", VariationType.SNP, "X",
                                          49111949, "G", "T",
                                          {"NM_001114377": MutationSyntax( \
                                              "NM_001114377", 756, 217, "", "")
                                          }, False, False)

        self.non_frame_shift_del = Variant("COSM1122495", VariationType.DEL, "X",
                                           49113232, "CTT", "",
                                           {"NM_001114377": MutationSyntax( \
                                               "NM_001114377", 615, 205, "", "")
                                           }, True, False)

        self.syn_homo_snp = Variant("COSM1122494", VariationType.SNP, "X",
                                    49112257, "C", "T",
                                    {"NM_001114377": MutationSyntax( \
                                        "NM_001114377", 653, 217, "", "")
                                    }, False, True)

        self.db_adapter = MartsAdapter(biomart="http://grch37.ensembl.org/biomart/martservice?query=")

    def test__incorp_snp(self):
        ts = list("TESTSEQUENCE")
        self.assertEqual(Generator._incorp_snp(ts, var_2, "tsc_1", 6, 6), 6)

    def test__incorp_insertion(self):
        ts = list("TESTSEQUENCE")
        self.assertEqual(Generator._incorp_insertion(ts, var_3, "tsc_1", 0, 0),  2)

    def test__incorp_deletion(self):
        ts = list("TESTSEQUEASDFGNCES")
        self.assertEqual(Generator._incorp_deletion(ts, var_4, "tsc_1", 0, 0),  -5)
        self.assertEqual(Generator._incorp_deletion(ts, var_6, "tsc_1", 0, 0),  -2)

    def test__check_for_problematic_variants(self):
        self.assertTrue(Generator._check_for_problematic_variants([var_2, var_1]))
        self.assertFalse(Generator._check_for_problematic_variants([var_5, var_6]))

    def test_non_syn_hetero_snp_trans_number(self):
        """
        tests if the number of generated transcripts for a heterozygous
        transcript is correct

        1 hetero vars = 2 transcripts
        :return:
        """
        vars_ = \
            [self.non_syn_hetero_snp, self.non_frame_shift_del, self.syn_homo_snp]

        trans = \
            [t for t in Generator.generate_transcripts_from_variants(vars_, self.db_adapter, EIdentifierTypes.REFSEQ)]

        self.assertTrue(len(trans) == 2**sum(not v.isHomozygous for v in vars_))

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
        dummy_vars = [var_3]
        trans = next(Generator.generate_transcripts_from_variants(dummy_vars, dummy_db, EIdentifierTypes.REFSEQ))
        self.assertEqual(str(trans), "AAAAACCTTCCCGGGGG")

        # SNPs:
        dummy_vars = [var_1]
        trans = next(Generator.generate_transcripts_from_variants(dummy_vars, dummy_db, EIdentifierTypes.REFSEQ))
        self.assertEqual(str(trans), "ATAAACCCCCGGGGG")

        # DELETIONS:
        dummy_vars = [var_4]
        trans = next(Generator.generate_transcripts_from_variants(dummy_vars, dummy_db, EIdentifierTypes.REFSEQ))
        self.assertEqual(str(trans), "AAAAAGGGGG")

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
        dummy_vars = [var_3, var_7, var_6]
        trans = next(Generator.generate_transcripts_from_variants(dummy_vars, dummy_db, EIdentifierTypes.REFSEQ))

        self.assertEqual(str(trans), "AAAAACCTTCTGGGG")

        # 2.) INS, DEL, INS
        dummy_vars = [var_9, var_4, var_8]
        trans = next(Generator.generate_transcripts_from_variants(dummy_vars, dummy_db, EIdentifierTypes.REFSEQ))
        self.assertEqual(str(trans), "AATTAAAGGGGGTTT")

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
        AAAAATTTCCGGGG (INS,DEL)
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
        trans_gener = Generator.generate_transcripts_from_variants(dummy_vars, dummy_db, EIdentifierTypes.REFSEQ)
        trans = [t for t in trans_gener]

        trans = list(map(str, trans))

        self.assertEqual(len(trans), 8)

        self.assertTrue("AAATTTGGGGG" in trans)
        self.assertTrue("AAAAATTTGGGGG" in trans)
        self.assertTrue("AAATTTCCCCCGGGGG" in trans)
        self.assertTrue("AAAAATTTCCCCCGGGGG" in trans)

        self.assertTrue("GGGTTTAAAAA" in trans)
        self.assertTrue("GGGGGTTTAAAAA" in trans)
        self.assertTrue("GGGTTTCCCCCAAAAA" in trans)
        self.assertTrue("GGGGGTTTCCCCCAAAAA" in trans)

        # def test_varinat_reader(self):
        #     vars = read_annovar_exonic("../IO/snp_annot_donor.txt.exonic_variant_function", gene_filter=[])
        #     print vars, len(vars)
        #     trans = list(generate_transcripts_from_variants(vars, self.db_adapter))
        #     print trans
        #     transToVar = {}
        #     for v in vars:
        #         for trans_id in v.coding.iterkeys():
        #             transToVar.setdefault(trans_id, []).append(v)

    #TODO: tests have to be reimplemented
    # def test_peptides_from_variants(self):
    #     """
    #     Create multiple peptides, given a set
    #     containing heterozygous variants .
    #
    #     Variants:
    #     3-DEL(-2)  , 5-INS(+3)  , 7-DEL(-4)
    #     HET-DEL(-2), HOM-INS(+3), HET-DEL(-1)
    #
    #     Reference sequence:
    #     AAAAACCCCCGGGGG
    #     AAATTTGGGGG (DEL,INS,DEL)
    #     AAATTTCCCCCGGGGG (DEL,INS)
    #     AAAAATTTGGGGG (INS,DEL)
    #     AAAAATTTCCCCCGGGGG (INS)
    #
    #     GGGGGCCCCCAAAAA
    #     GGGTTTCAAAAA (DEL,INS,DEL)
    #     GGGTTTCCCCCAAAAA (DEL,INS)
    #     GGGGGTTTCAAAAA (INS,DEL)
    #     GGGGGTTTCCCCCAAAAA (INS)
    #
    #
    #     Resulting protein sequences:
    #     KFG
    #     KNLG
    #     KFPPG
    #     KNFPRG
    #
    #     GFK
    #     GGLK
    #     GFPPK
    #     GGFPQK
    #
    #     Resulting peptides of length 3:
    #     KFG +
    #     KNL +
    #     NLG +
    #     KFP +
    #     FPP +
    #     PPG +
    #     KNF +
    #     NFP +
    #     FPR +
    #     PRG +
    #
    #     GFK +
    #     GGL +
    #     GLK +
    #     GFP +
    #     FPP +
    #     PPK +
    #     GGF +
    #     GFP +
    #     FPQ +
    #     PQK +
    #     """
    #     dummy_db = DummyAdapter()
    #
    #     exp_peps = set(['PRG', 'GLK', 'PPG', 'KFP', 'GFK', 'PPK', 'GFP', 'PQK', 'KNL', 'KFG', 'GGF', 'FPQ',
    #                     'FPP', 'NLG', 'FPR', 'KNF', 'GGL', 'NFP'])
    #     # 1) INS, SNP, DEL
    #     dummy_vars = [var_10, var_11, var_12]
    #     peps = set(map(lambda x: str(x), Generator.generate_peptides_from_variants(dummy_vars, 3, dummy_db, EIdentifierTypes.REFSEQ)))
    #
    #     peps_from_prot = set(map(str, Generator.generate_peptides_from_proteins(Generator.generate_proteins_from_transcripts(
    #                        Generator.generate_transcripts_from_variants(dummy_vars, dummy_db, EIdentifierTypes.REFSEQ)),
    #         3)))
    #
    #     self.assertTrue(len(peps - peps_from_prot) == 0)
    #     self.assertTrue(len(peps_from_prot - peps) == 0)
    #     self.assertTrue(len(peps-exp_peps) == 0)
    #     self.assertTrue(len(exp_peps-peps) == 0)

    # def test_real_life_test(self):
    #     mart = MartsAdapter(biomart="http://grch37.ensembl.org/biomart/martservice?query=")
    #
    #     ano_path = os.path.join(os.path.dirname(inspect.getfile(Fred2)), "Data/examples/test_annovar.out")
    #     vars = read_annovar_exonic(ano_path)
    #
    #     peps = set(map(lambda x: str(x), Generator.generate_peptides_from_variants(vars, 9, mart,
    #                                                                                EIdentifierTypes.REFSEQ)))
    #
    #     peps_from_prot = set(map(str, Generator.generate_peptides_from_proteins(
    #         Generator.generate_proteins_from_transcripts(
    #         Generator.generate_transcripts_from_variants(vars, mart, EIdentifierTypes.REFSEQ)), 9)))
    #
    #     self.assertTrue(len(peps - peps_from_prot) == 0)
    #     self.assertTrue(len(peps_from_prot - peps) == 0)

    def test_peptides_from_varaints(self):
        coding = {}
        coding['NM_080751'] = MutationSyntax('NM_080751',2629,876,'c.2630C>T','p.Pro877Leu')
        var = Variant('line0',0,20,2621905,'C','T',coding,True,False)
        var.gene = 'TMC2'
        ma = MartsAdapter(biomart="http://ensembl.org")

        vars = [var, Variant("testInsertion", 2, 20, 2621899, "", "AAAAAA", {'NM_080751':MutationSyntax('NM_080751',2625,876,'c.2630C>T','p.Pro877Leu')}, True, False)]

        test = Generator.generate_peptides_from_variants(vars, 9, ma, id_type=EIdentifierTypes.REFSEQ, peptides=None)
        test2 = [x for x in test]
        print((len(test2)))

        ts = list()
        #using a tweaked generator that takes another sequence source if the sequence is too short in respect to the given variants
        #in this case a newer/older sequence from mart in respect to what was given as reference in the annotation process
        t = Generator.generate_transcripts_from_variants(vars, ma, id_type=EIdentifierTypes.REFSEQ)
        ts = [x for x in t]
        print((len(ts[0])))
        p = Generator.generate_proteins_from_transcripts(ts, to_stop=True)
        ps = [x for x in p]
        e = Generator.generate_peptides_from_proteins(ps, 9)
        es = [x for x in e]
        print((len(es)))

        #print vars
        print(len(vars))


    def test_proteins_from_variants(self):
        """
                Variants:
        3-DEL(-2)  , 5-INS(+3)  , 7-DEL(-4)
        HET-DEL(-2), HOM-INS(+3), HET-DEL(-1)

        Reference sequence:set(['GLK', 'PPK', 'GFP', 'PQK', 'GFK', 'GGF', 'FPQ', 'FPP', 'GGL'])
        AAATTTCGGGGG (DEL,INS,DEL)
        AAATTTCCCCCGGGGG (DEL,INS)
        AAAAATTTCCGGGG (INS,DEL)
        AAAAATTTCCCCCGGGGG (INS)

        GGGTTTCAAAAA (DEL,INS,DEL)
        GGGTTTCCCCCAAAAA (DEL,INS)
        GGGGGTTTCAAAAA (INS,DEL)
        GGGGGTTTCCCCCAAAAA (INS)


        Resulting protein sequences:
        KFG
        KNLG
        KFPPG
        KNFPRG

        GFK
        GGLK
        GFPPK
        GGFPQK
        """
        dummy_db = DummyAdapter()
        dummy_vars = [var_10, var_11, var_12]
        exp_prot = set(['KFG', 'KNLG', 'KFPPG', 'KNFPRG', 'GFK', 'GGLK', 'GFPPK', 'GFPPK', 'GGFPQK'])
        prot = set([str(x) for x in Generator.generate_proteins_from_transcripts(
                           Generator.generate_transcripts_from_variants(dummy_vars, dummy_db, EIdentifierTypes.REFSEQ))]
                   )
        self.assertTrue(len(prot-exp_prot) == 0)
        self.assertTrue(len(exp_prot-prot) == 0)
