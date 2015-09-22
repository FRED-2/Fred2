from unittest import TestCase

__author__ = 'walzer,schubert'

from Fred2.Core import Transcript
from Fred2.Core import Variant
from Fred2.Core import VariationType
from Fred2.Core import MutationSyntax
from Fred2.Core import generate_proteins_from_transcripts


class TestTranscript(TestCase):
    def setUp(self):
        self.simple = Transcript("")
        self.simple_new = Transcript("")
        self.w_gid = Transcript("", _gene_id="123")
        self.w_tid = Transcript("", _transcript_id="tid")
        self.w_id = Transcript("", "gid", "tid")

        #Internal indexing starts at 0! MutationSyntax coming e.g. from ANNOVAR starts at 1!
        self.gcg_v1 = Variant("rs5650", VariationType.SNP, 2, 162145588, 'G', 'T',
                         {"NM_002054.4": MutationSyntax("NM_002054.4", 343, 114, "c.344C>A", "p.A115D")}, False, False)

        self.gcg_v1.offsets["NM_002054.4"] = 0
        self.gcg_ts = "gcatagaatgcagatgagcaaagtgagtgggagagggaagtcatttgtaacaaaaactcattatttacagatgagaaatttatattgtcagcgtaatatctgtgaggctaaacagagctggagagtatataaaagcagtgcgccttggtgcagaagtacagagcttaggacacagagcacatcaaaagttcccaaagagggcttgctctctcttcacctgctctgttctacagcacactaccagaagacagcagaaatgaaaagcatttactttgtggctggattatttgtaatgctggtacaaggcagctggcaacgttcccttcaagacacagaggagaaatccagatcattctcagcttcccaggcagacccactcagtgatcctgatcagatgaacgaggacaagcgccattcacagggcacattcaccagtgactacagcaagtatctggactccaggcgtgcccaagattttgtgcagtggttgatgaataccaagaggaacaggaataacattgccaaacgtcacgatgaatttgagagacatgctgaagggacctttaccagtgatgtaagttcttatttggaaggccaagctgccaaggaattcattgcttggctggtgaaaggccgaggaaggcgagatttcccagaagaggtcgccattgttgaagaacttggccgcagacatgctgatggttctttctctgatgagatgaacaccattcttgataatcttgccgccagggactttataaactggttgattcagaccaaaatcactgacaggaaataactatatcactattcaagatcatcttcacaacatcacctgctagccacgtgggatgtttgaaatgttaagtcctgtaaatttaagaggtgtattctgaggccacattgctttgcatgccaataaataaattttcttttagtgttgtgtagccaaaaattacaaatggaataaagttttatcaaaatattgctaaaatatcagctttaaaatatgaaagtgctagattctgttattttcttcttattttggatgaagtaccccaacctgtttacatttagcgataaaattatttttctatgatataatttgtaaatgtaaattattccgatctgacatatctgcattataataataggagaatagaagaactggtagccacagtggtgaaattggaaagagaactttcttcctgaaacctttgtcttaaaaatactcagctttcaatgtatcaaagatacaattaaataaaattttcaagcttctttaccattgtct"
        self.w_v = Transcript(self.gcg_ts, 'GLUC_HUMAN', "NM_002054.4", [self.gcg_v1])

    def test_consistency(self):
        self.assertTrue(repr(self.simple) == "TRANSCRIPT: 0\n\tVARIANTS:\n\tSEQUENCE:  (mRNA)")
        self.assertTrue(repr(self.simple_new) == "TRANSCRIPT: 1\n\tVARIANTS:\n\tSEQUENCE:  (mRNA)")
        self.assertTrue(repr(self.w_gid) == "TRANSCRIPT: 2\n\tVARIANTS:\n\tSEQUENCE:  (mRNA)")
        self.assertTrue(repr(self.w_tid) == "TRANSCRIPT: tid\n\tVARIANTS:\n\tSEQUENCE:  (mRNA)")
        self.assertTrue(repr(self.w_id) == "TRANSCRIPT: tid\n\tVARIANTS:\n\tSEQUENCE:  (mRNA)")
        self.assertTrue(repr(
            self.w_v) == "TRANSCRIPT: NM_002054.4\n\tVARIANTS:\n\t\tpos 343: Variant(g.162145588G>T)\n\tSEQUENCE: GCATAGAATGCAGATGAGCAAAGTGAGTGGGAGAGGGAAGTCATTTGTAACAAAAACTCATTATTTACAGATGAGAAATTTATATTGTCAGCGTAATATCTGTGAGGCTAAACAGAGCTGGAGAGTATATAAAAGCAGTGCGCCTTGGTGCAGAAGTACAGAGCTTAGGACACAGAGCACATCAAAAGTTCCCAAAGAGGGCTTGCTCTCTCTTCACCTGCTCTGTTCTACAGCACACTACCAGAAGACAGCAGAAATGAAAAGCATTTACTTTGTGGCTGGATTATTTGTAATGCTGGTACAAGGCAGCTGGCAACGTTCCCTTCAAGACACAGAGGAGAAATCCAGATCATTCTCAGCTTCCCAGGCAGACCCACTCAGTGATCCTGATCAGATGAACGAGGACAAGCGCCATTCACAGGGCACATTCACCAGTGACTACAGCAAGTATCTGGACTCCAGGCGTGCCCAAGATTTTGTGCAGTGGTTGATGAATACCAAGAGGAACAGGAATAACATTGCCAAACGTCACGATGAATTTGAGAGACATGCTGAAGGGACCTTTACCAGTGATGTAAGTTCTTATTTGGAAGGCCAAGCTGCCAAGGAATTCATTGCTTGGCTGGTGAAAGGCCGAGGAAGGCGAGATTTCCCAGAAGAGGTCGCCATTGTTGAAGAACTTGGCCGCAGACATGCTGATGGTTCTTTCTCTGATGAGATGAACACCATTCTTGATAATCTTGCCGCCAGGGACTTTATAAACTGGTTGATTCAGACCAAAATCACTGACAGGAAATAACTATATCACTATTCAAGATCATCTTCACAACATCACCTGCTAGCCACGTGGGATGTTTGAAATGTTAAGTCCTGTAAATTTAAGAGGTGTATTCTGAGGCCACATTGCTTTGCATGCCAATAAATAAATTTTCTTTTAGTGTTGTGTAGCCAAAAATTACAAATGGAATAAAGTTTTATCAAAATATTGCTAAAATATCAGCTTTAAAATATGAAAGTGCTAGATTCTGTTATTTTCTTCTTATTTTGGATGAAGTACCCCAACCTGTTTACATTTAGCGATAAAATTATTTTTCTATGATATAATTTGTAAATGTAAATTATTCCGATCTGACATATCTGCATTATAATAATAGGAGAATAGAAGAACTGGTAGCCACAGTGGTGAAATTGGAAAGAGAACTTTCTTCCTGAAACCTTTGTCTTAAAAATACTCAGCTTTCAATGTATCAAAGATACAATTAAATAAAATTTTCAAGCTTCTTTACCATTGTCT (mRNA)")

    def test_translate(self):
        gcg_var = "A*NADEQSEWEREVICNKNSLFTDEKFILSA*YL*G*TELESI*KQCALVQKYRA*DTEHIKSSQRGLALSSPALFYSTLPEDSRNEKHLLCGWIICNAGTRQLATFPSRHRGEIQIILSFPGRPTQ*S*SDERGQAPFTGHIHQ*LQQVSGLQACPRFCAVVDEYQEEQE*HCQTSR*I*ETC*RDLYQ*CKFLFGRPSCQGIHCLAGERPRKARFPRRGRHC*RTWPQTC*WFFL**DEHHS**SCRQGLYKLVDSDQNH*QEITISLFKIIFTTSPASHVGCLKC*VL*I*EVYSEATLLCMPINKFSFSVV*PKITNGIKFYQNIAKISALKYESARFCYFLLILDEVPQPVYI*R*NYFSMI*FVNVNYSDLTYLHYNNRRIEELVATVVKLERELSS*NLCLKNTQLSMYQRYN*IKFSSFFTIV"
        self.assertTrue(str(generate_proteins_from_transcripts(self.w_v).next()) == gcg_var)
        # http://stackoverflow.com/questions/3892218/how-to-test-with-pythons-unittest-that-a-warning-has-been-thrown

    def test_indexing(self):
        self.assertTrue(self.w_v[0] == "G")

    def test_transcript_slicing(self):
        new_t = self.w_v[343:]
        self.assertTrue(self.gcg_v1.get_transcript_position(new_t.transcript_id) == 0)

