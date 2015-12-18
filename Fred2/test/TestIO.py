from unittest import TestCase
import copy

from Fred2.Core import Allele
from Fred2.IO import FileReader
from Fred2.IO.MartsAdapter import MartsAdapter
from Fred2.IO.EnsemblAdapter import EnsemblDB
from Fred2.IO.RefSeqAdapter import RefSeqAdapter
from Fred2.IO.UniProtAdapter import UniProtDB
from Fred2.IO.ADBAdapter import ADBAdapter, EAdapterFields, EIdentifierTypes
import warnings
import logging
import os
import inspect
import Fred2

__author__ = 'walzer'


class TestIO(TestCase):
    def assertWarnings(self, warning, call, *args, **kwds):
        with warnings.catch_warnings(record=True) as warning_list:
            warnings.simplefilter('always')
            result = call(*args, **kwds)
            self.assertTrue(any(item.category == warning for item in warning_list))

    def setUp(self):
        self.ale_path = os.path.join(os.path.dirname(inspect.getfile(Fred2)), "Data/examples/alleles.txt")
        self.ale_zonk_path = os.path.join(os.path.dirname(inspect.getfile(Fred2)), "Data/examples/alleles_defect.txt")
        self.ale_no_path = os.path.join(os.path.dirname(inspect.getfile(Fred2)), "Data/examples/magic.txt")
        self.fa_path = os.path.join(os.path.dirname(inspect.getfile(Fred2)), "Data/examples/testSequences.fasta")
        self.fa_unconventional_path = os.path.join(os.path.dirname(inspect.getfile(Fred2)), "Data/examples/testSequences_d2s.fasta")
        self.edb_cds_path = os.path.join(os.path.dirname(inspect.getfile(Fred2)), "Data/examples/Homo_sapiens.GRCh38.cds.test_stub.fa")
        self.edb_pep_path = os.path.join(os.path.dirname(inspect.getfile(Fred2)), "Data/examples/Homo_sapiens.GRCh38.pep.test_stub.fa")
        self.ano_path = os.path.join(os.path.dirname(inspect.getfile(Fred2)), "Data/examples/test_annovar.out")
        self.NP_001005353 = "MASKLLRAVILGPPGSGKGTVCQRIAQNFGLQHLSSGHFLRENIKASTEVGEMAKQYIEKSLLVPDHVITRLMMSELENRRGQHWLLDGFPRTLGQAEALDKICEVDLVISLNIPFETLKDRLSRRWIHPPSGRVYNLDFNPPHVHGIDDVTGEPLVQQEDDKPEAVAARLRQYKDVAKPVIELYKSRGVLHQFSGTETNKIWPYVYTLFSNKITPIQSKEAY"
        self.ENSP00000369497 = "MPIGSKERPTFFEIFKTRCNKADLGPISLNWFEELSSEAPPYNSEPAEESEHKNNNYEPNLFKTPQRKPSYNQLASTPIIFKEQGLTLPLYQSPVKELDKFKLDLGRNVPNSRHKSLRTVKTKMDQADDVSCPLLNSCLSESPVVLQCTHVTPQRDKSVVCGSLFHTPKFVKGRQTPKHISESLGAEVDPDMSWSSSLATPPTLSSTVLIVRNEEASETVFPHDTTANVKSYFSNHDESLKKNDRFIASVTDSENTNQREAASHGFGKTSGNSFKVNSCKDHIGKSMPNVLEDEVYETVVDTSEEDSFSLCFSKCRTKNLQKVRTSKTRKKIFHEANADECEKSKNQVKEKYSFVSEVEPNDTDPLDSNVANQKPFESGSDKISKEVVPSLACEWSQLTLSGLNGAQMEKIPLLHISSCDQNISEKDLLDTENKRKKDFLTSENSLPRISSLPKSEKPLNEETVVNKRDEEQHLESHTDCILAVKQAISGTSPVASSFQGIKKSIFRIRESPKETFNASFSGHMTDPNFKKETEASESGLEIHTVCSQKEDSLCPNLIDNGSWPATTTQNSVALKNAGLISTLKKKTNKFIYAIHDETSYKGKKIPKDQKSELINCSAQFEANAFEAPLTFANADSGLLHSSVKRSCSQNDSEEPTLSLTSSFGTILRKCSRNETCSNNTVISQDLDYKEAKCNKEKLQLFITPEADSLSCLQEGQCENDPKSKKVSDIKEEVLAAACHPVQHSKVEYSDTDFQSQKSLLYDHENASTLILTPTSKDVLSNLVMISRGKESYKMSDKLKGNNYESDVELTKNIPMEKNQDVCALNENYKNVELLPPEKYMRVASPSRKVQFNQNTNLRVIQKNQEETTSISKITVNPDSEELFSDNENNFVFQVANERNNLALGNTKELHETDLTCVNEPIFKNSTMVLYGDTGDKQATQVSIKKDLVYVLAEENKNSVKQHIKMTLGQDLKSDISLNIDKIPEKNNDYMNKWAGLLGPISNHSFGGSFRTASNKEIKLSEHNIKKSKMFFKDIEEQYPTSLACVEIVNTLALDNQKKLSKPQSINTVSAHLQSSVVVSDCKNSHITPQMLFSKQDFNSNHNLTPSQKAEITELSTILEESGSQFEFTQFRKPSYILQKSTFEVPENQMTILKTTSEECRDADLHVIMNAPSIGQVDSSKQFEGTVEIKRKFAGLLKNDCNKSASGYLTDENEVGFRGFYSAHGTKLNVSTEALQKAVKLFSDIENISEETSAEVHPISLSSSKCHDSVVSMFKIENHNDKTVSEKNNKCQLILQNNIEMTTGTFVEEITENYKRNTENEDNKYTAASRNSHNLEFDGSDSSKNDTVCIHKDETDLLFTDQHNICLKLSGQFMKEGNTQIKEDLSDLTFLEVAKAQEACHGNTSNKEQLTATKTEQNIKDFETSDTFFQTASGKNISVAKESFNKIVNFFDQKPEELHNFSLNSELHSDIRKNKMDILSYEETDIVKHKILKESVPVGTGNQLVTFQGQPERDEKIKEPTLLGFHTASGKKVKIAKESLDKVKNLFDEKEQGTSEITSFSHQWAKTLKYREACKDLELACETIEITAAPKCKEMQNSLNNDKNLVSIETVVPPKLLSDNLCRQTENLKTSKSIFLKVKVHENVEKETAKSPATCYTNQSPYSVIENSALAFYTSCSRKTSVSQTSLLEAKKWLREGIFDGQPERINTADYVGNYLYENNSNSTIAENDKNHLSEKQDTYLSNSSMSNSYSYHSDEVYNDSGYLSKNKLDSGIEPVLKNVEDQKNTSFSKVISNVKDANAYPQTVNEDICVEELVTSSSPCKNKNAAIKLSISNSNNFEVGPPAFRIASGKIVCVSHETIKKVKDIFTDSFSKVIKENNENKSKICQTKIMAGCYEALDDSEDILHNSLDNDECSTHSHKVFADIQSEEILQHNQNMSGLEKVSKISPCDVSLETSDICKCSIGKLHKSVSSANTCGIFSTASGKSVQVSDASLQNARQVFSEIEDSTKQVFSKVLFKSNEHSDQLTREENTAIRTPEHLISQKGFSYNVVNSSAFSGFSTASGKQVSILESSLHKVKGVLEEFDLIRTEHSLHYSPTSRQNVSKILPRVDKRNPEHCVNSEMEKTCSKEFKLSNNLNVEGGSSENNHSIKVSPYLSQFQQDKQQLVLGTKVSLVENIHVLGKEQASPKNVKMEIGKTETFSDVPVKTNIEVCSTYSKDSENYFETEAVEIAKAFMEDDELTDSKLPSHATHSLFTCPENEEMVLSNSRIGKRRGEPLILVGEPSIKRNLLNEFDRIIENQEKSLKASKSTPDGTIKDRRLFMHHVSLEPITCVPFRTTKERQEIQNPNFTAPGQEFLSKSHLYEHLTLEKSSSNLAVSGHPFYQVSATRNEKMRHLITTGRPTKVFVPPFKTKSHFHRVEQCVRNINLEENRQKQNIDGHGSDDSKNKINDNEIHQFNKNNSNQAVAVTFTKCEEEPLDLITSLQNARDIQDMRIKKKQRQRVFPQPGSLYLAKTSTLPRISLKAAVGGQVPSACSHKQLYTYGVSKHCIKINSKNAESFQFHTEDYFGKESLWTGKGIQLADGGWLIPSNDGKAGKEEFYRALCDTPGVDPKLISRIWVYNHYRWIIWKLAAMECAFPKEFANRCLSPERVLLQLKYRYDTEIDRSRRSAIKKIMERDDTAAKTLVLCVSDIISLSANISETSSNKTSSADTQKVAIIELTDGWYAVKAQLDPPLLAVLKNGRLTVGQKIILHGAELVGSPDACTPLEAPESLMLKISANSTRPARWYTKLGFFPDPRPFPLPLSSLFSDGGNVGCVDVIIQRAYPIQWMEKTSSGLYIFRNEREEEKEAAKYVEAQQKRLEALFTKIQEEFEEHEENTTKPYLPSRALTRQQVRALQDGAELYEAVKNAADPAYLEGYFSEEQLRALNNHRQMLNDKKQAQIQLEIRKAMESAEQKEQGLSRDVTTVWKLRIVSYSKKEKDSVILSIWRPSSDLYSLLTEGKRYRIYHLATSKSKSKSERANIQLAATKKTQYQQLPVSDEILFQIYQPREPLHFSKFLDPDFQPSCSEVDLIGFVVSVVKKTGLAPFVYLSDECYNLLAIKFWIDLNEDIIKPHMLIAASNLQWRPESKSGLLTLFAGDFSVFSASPKEGHFQETFNKMKNTVENIDILCNEAENKLMHILHANDPKWSTPTKDCTSGPYTAQIIPGTGNKLLMSSPNCEIYYQSPLSLCMAKRKSVSTPVSAQMTSKSCKGEKEIDDQKNCKKRRALDFLSRLPLPPPVSPICTFVSPAAQKAFQPPRSCGTKYETPIKKKELNSPQMTPFKKFNEISLLESNSIADEELALINTQALLSGSTGEKQFISVSESTRTAPTSSEDYLRLKRRCTTSLIKEQESSQASTEECEKNKQDTITTKKYI"
        self.ENST00000361221 = {0: '', 1: '+', 2: 'ATGGCTTCCTCCAACCCTCCTCCACAGCCTGCCATAGGAGATCAGCTGGTTCCAGGAGTCCCAGGCCCCTCCTCTGAGGCAGAGGACGACCCAGGAGAGGCGTTTGAGTTTGATGACAGTGATGATGAAGAGGACACCAGCGCAGCCCTGGGCGTCCCCAGCCTTGCTCCTGAGAGGGACACAGACCCCCCACTGATCCACTTGGACTCCATCCCTGTCACTGACCCAGACCCAGCAGCTGCTCCACCCGGCACAGGGGTGCCAGCCTGGGTGAGCAATGGGGATGCAGCGGACGCAGCCTTCTCCGGGGCCCGGCACTCCAGCTGGAAGCGGAAGAGTTCCCGTCGCATTGACCGGTTCACTTTCCCCGCCCTGGAAGAGGATGTGATTTATGACGACGTCCCCTGCGAGAGCCCAGATGCGCATCAGCCCGGGGCAGAGAGGAACCTGCTCTACGAGGATGCGCACCGGGCTGGGGCCCCTCGGCAGGCGGAGGACCTAGGCTGGAGCTCCAGTGAGTTCGAGAGCTACAGCGAGGACTCGGGGGAGGAGGCCAAGCCGGAGGTCGAGGTCGAGCCCGCCAAGCACCGAGTGTCCTTCCAGCCCAAGCTTTCTCCAGACCTGACTAGGCTAAAGGAGAGATACGCCAGGACTAAGAGAGACATCTTGGCTTTGAGAGTTGGGGGGAGAGACATGCAGGAGCTGAAGCACAAGTACGATTGTAAGATGACCCAGCTCATGAAGGCCGCCAAGAGCGGGACCAAGGATGGGCTGGAGAAGACACGGATGGCCGTGATGCGCAAAGTCTCCTTCCTGCACAGGAAGGACGTCCTCGGTGACTCGGAGGAGGAGGACATGGGGCTCCTGGAGGTCAGCGTTTCGGACATCAAGCCCCCAGCCCCAGAGCTGGGCCCCATGCCAGAGGGCCTGAGCCCTCAGCAGGTGGTCCGGAGGCATATCCTGGGCTCCATCGTGCAGAGCGAAGGCAGCTACGTGGAGTCTCTGAAGCGGATACTCCAGGACTACCGCAACCCCCTGATGGAGATGGAGCCCAAGGCGCTGAGCGCCCGCAAGTGCCAGGTGGTGTTCTTCCGCGTGAAGGAGATCCTGCACTGCCACTCCATGTTCCAGATCGCCCTGTCCTCCCGCGTGGCTGAGTGGGATTCCACCGAGAAGATCGGGGACCTCTTCGTGGCCTCGTTTTCCAAGTCCATGGTGCTAGATGTGTACAGTGACTACGTGAACAACTTCACCAGTGCCATGTCCATCATCAAGAAGGCCTGCCTCACCAAGCCTGCCTTCCTCGAGTTCCTCAAGCGACGGCAGGTGTGCAGCCCAGACCGTGTCACCCTCTACGGGCTGATGGTCAAGCCCATCCAGAGGTTCCCACAGTTCATACTCCTGCTTCAGGACATGCTGAAGAACACCCCCAGGGGCCATCCGGACAGGCTGTCGCTGCAGCTGGCCCTCACAGAGCTGGAGACGCTGGCTGAGAAGCTGAACGAGCAGAAGCGGCTGGCTGACCAGGTGGCTGAGATCCAGCAGCTGACCAAGAGCGTCAGTGACCGCAGCAGCCTCAACAAGCTGTTGACCTCAGGCCAGCGGCAGCTGCTCCTGTGTGAGACGTTGACGGAGACCGTGTACGGTGACCGCGGGCAGCTAATTAAGTCCAAGGAGCGTCGGGTCTTCCTGCTCAACGACATGCTTGTCTGTGCCAACATCAACTTCAAGCCTGCCAACCACAGGGGCCAGCTGGAGATCAGCAGCCTGGTGCCCCTGGGGCCCAAGTATGTGGTGAAGTGGAACACGGCGCTGCCCCAGGTGCAGGTGGTGGAGGTGGGCCAGGACGGTGGCACCTATGACAAGGACAATGTGCTCATCCAGCACTCAGGCGCCAAGAAGGCCTCTGCCTCAGGGCAGGCTCAGAATAAGGTGTACCTCGGCCCCCCACGCCTCTTCCAGGAGCTGCAGGACCTGCAGAAGGACCTGGCCGTGGTGGAGCAGATCACGCTTCTCATCAGCACGCTGCACGGCACCTACCAGAACCTGAACATGACTGTGGCTCAAGACTGGTGCCTGGCCCTGCAGAGGCTGATGCGGGTGAAGGAGGAAGAGATCCACTCGGCCAACAAGTGCCGTCTCAGGCTCCTGCTTCCTGGGAAACCCGACAAGTCCGGCCGCCCCATTAGCTTCATGGTGGTTTTCATCACCCCCAACCCCCTGAGCAAGATTTCCTGGGTCAACAGGTTACATTTGGCCAAAATCGGACTCCGGGAGGAGAACCAGCCAGGCTGGCTATGCCCGGATGAGGACAAGAAGAGCAAAGCCCCATTCTGGTGCCCGATCCTGGCCTGCTGCATCCCTGCCTTCTCCTCCCGGGCACTCAGCCTGCAGCTTGGGGCCCTGGTCCACAGTCCTGTCAACTGTCCCCTGCTGGGTTTCTCAGCAGTCAGCACCTCCCTTCCACAGGGCTACCTCTGGGTCGGGGGCGGACAGGAAGGCGCAGGGGGCCAGGTGGAAATCTTTTCCTTGAACCGGCCCTCGCCCCGCACCGTCAAGTCCTTCCCACTGGCAGCCCCTGTGCTCTGCATGGAGTATATCCCGGAGCTGGAGGAGGAGGCGGAGAGCAGAGACGAGAGCCCGACAGTTGCTGACCCCTCGGCCACGGTGCATCCAACCATCTGCCTCGGGCTCCAGGATGGCAGCATCCTCCTCTACAGCAGTGTGGACACTGGCACCCAGTGCCTGGTGAGCTGCAGGAGCCCAGGTCTGCAGCCTGTGCTCTGCCTGCGACACAGCCCCTTCCACCTGCTCGCTGGCCTGCAGGATGGGACCCTTGCTGCTTACCCTCGGACCAGCGGAGGTGTCCTGTGGGACCTGGAGAGCCCTCCCGTGTGCCTGACTGTGGGGCCCGGGCCTGTCCGCACCCTGTTGAGCCTGGAGGATGCCGTGTGGGCCAGCTGTGGGCCCTGGGTCACTGTCCTGGAAGCCACCACCCTGCAGCCTCAGCAAAGCTTCGAGGCGCACCAGGACGAGGCAGTGAGCGTGACACACATGGTGAAGGCGGGCAGCGGCGTCTGGATGGCCTTCTCCTCCGGCACCTCCATCCGCCTCTTCCACACTGAGACCCTGGAGCATCTGCAAGAGATCAACATCGCCACCAGGACCACCTTCCTCCTGCCAGGCCAGAAGCACTTGTGTGTCACCAGCCTCCTGATCTGCCAGGGTCTGCTCTGGGTGGGCACTGACCAGGGTGTCATCGTCCTGCTGCCCGTGCCTCGGCTGGAAGGCATCCCCAAGATCACAGGGAAAGGCATGGTCTCACTCAACGGGCACTGTGGGCCTGTGGCCTTCCTGGCTGTGGCTACCAGCATCCTGGCCCCTGACATCCTGCGGAGTGACCAGGAGGAGGCTGAGGGGCCCCGGGCTGAGGAGGACAAGCCAGACGGGCAGGCACACGAGCCCATGCCCGATAGCCACGTGGGCCGAGAGCTGACCCGCAAGAAGGGCATCCTCTTGCAGTACCGCCTGCGCTCCACCGCACACCTCCCGGGCCCGCTGCTCTCCATGCGGGAGCCGGCGCCTGCTGATGGCGCAGCTTTGGAGCACAGCGAGGAGGACGGCTCCATTTACGAGATGGCCGACGACCCCGACATCTGGGTGCGCAGCCGGCCCTGCGCCCGCGACGCCCACCGCAAGGAGATTTGCTCTGTGGCCATCATCTCCGGCGGGCAGGGCTACCGCAACTTTGGCAGCGCTCTGGGCAGCAGTGGGAGGCAGGCCCCGTGTGGGGAGACGGACAGCACCCTCCTCATCTGGCAGGTGCCCTTGATGCTATAG'}
        self.ENSEMBL_ensg = "ENSG00000138674"

    def test_read_lines(self):
        alleles = FileReader.read_lines(self.ale_path, in_type=Allele)
        self.assertEqual(len(alleles), 2)
        self.assertRaises(IOError, FileReader.read_lines, self.ale_no_path, in_type=Allele)
        self.assertRaises(ValueError, FileReader.read_lines, self.ale_zonk_path, in_type=Allele)

    def test_read_fasta(self):
        seqs = FileReader.read_fasta(self.fa_path)
        self.assertEqual(len(seqs), 2)
        self.assertRaises(IndexError, FileReader.read_fasta, self.fa_unconventional_path) # no "|"

    def test_read_annovar_exonic(self):
        ano = FileReader.read_annovar_exonic(self.ano_path)
        self.assertEqual(len(ano), 5)

    def test_EnsemblAdapter(self):
        ed = EnsemblDB()
        ed.read_seqs(self.edb_cds_path)
        self.assertEqual(len(ed.collection), 30)
        self.assertEqual(ed.get_transcript_sequence("ENST00000348405", type=EIdentifierTypes.ENSEMBL),
                         ed.get_transcript_information("ENST00000348405", type=EIdentifierTypes.ENSEMBL)[2])

        ed.read_seqs(self.edb_pep_path)
        self.assertEqual(len(ed.collection), 60)
        self.assertEqual(ed.get_product_sequence("ENSP00000337602", type=EIdentifierTypes.ENSEMBL).seq,
                         ed.get_transcript_information("ENSP00000337602", type=EIdentifierTypes.ENSEMBL)[2])
        self.assertEqual(ed.get_transcript_information("ENSP00000337602", type=EIdentifierTypes.ENSEMBL)[1], '-')
        self.assertEqual(ed.get_transcript_information("ENSP00000337602", type=EIdentifierTypes.ENSEMBL)[0],
                         self.ENSEMBL_ensg)

    def test_MartsAdapter(self):
        ma = MartsAdapter(biomart="http://grch37.ensembl.org")

        self.assertEqual((1515, 1529), ma.get_transcript_position('ENST00000361221', '17953929', '17953943', type=EIdentifierTypes.ENSEMBL))
        self.assertIsNone(ma.get_transcript_position("ENST00000614237", 7566927, 7566927, type=EIdentifierTypes.ENSEMBL))
        #logging.captureWarnings(True)
        #result = ma.get_transcript_position("ENST00000614237", 7566927, 7566927, type=EIdentifierTypes.ENSEMBL)
        self.assertEqual("TP53", ma.get_gene_by_position(17, 7566927, 7566927))
        self.assertIsNone(ma.get_product_sequence("Q15942", type=EIdentifierTypes.UNIPROT))
        self.assertEqual(self.NP_001005353, ma.get_product_sequence("NP_001005353", type=EIdentifierTypes.REFSEQ))
        self.assertEqual(self.ENSP00000369497, ma.get_product_sequence("ENSP00000369497", type=EIdentifierTypes.ENSEMBL))
        self.assertEqual(self.ENST00000361221[2], ma.get_transcript_sequence('ENST00000361221', type=EIdentifierTypes.ENSEMBL))
        self.assertIsNone(ma.get_transcript_sequence("ENST00000614237", type=EIdentifierTypes.ENSEMBL))
        self.assertDictEqual(self.ENST00000361221, ma.get_transcript_information('ENST00000361221', type=EIdentifierTypes.ENSEMBL))
        self.assertIsNone(ma.get_transcript_information("ENST00000614237", type=EIdentifierTypes.ENSEMBL))
        self.assertEqual(str(ma.get_ensembl_ids_from_id('TP53', type=EIdentifierTypes.GENENAME)), "[{0: 'ENSG00000141510', 1: '-', 3: 'ENST00000413465', 4: 'ENSP00000410739'}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000359597', 4: 'ENSP00000352610'}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000504290', 4: ''}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000510385', 4: ''}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000504937', 4: ''}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000269305', 4: 'ENSP00000269305'}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000455263', 4: 'ENSP00000398846'}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000420246', 4: 'ENSP00000391127'}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000445888', 4: 'ENSP00000391478'}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000576024', 4: 'ENSP00000458393'}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000509690', 4: 'ENSP00000425104'}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000514944', 4: 'ENSP00000423862'}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000574684', 4: ''}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000505014', 4: ''}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000508793', 4: 'ENSP00000424104'}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000604348', 4: 'ENSP00000473895'}, {0: 'ENSG00000141510', 1: '-', 3: 'ENST00000503591', 4: 'ENSP00000426252'}, {0: 'LRG_321', 1: '+', 3: 'LRG_321t8', 4: 'LRG_321p8'}, {0: 'LRG_321', 1: '+', 3: 'LRG_321t7', 4: 'LRG_321p13'}, {0: 'LRG_321', 1: '+', 3: 'LRG_321t6', 4: 'LRG_321p12'}, {0: 'LRG_321', 1: '+', 3: 'LRG_321t5', 4: 'LRG_321p11'}, {0: 'LRG_321', 1: '+', 3: 'LRG_321t4', 4: 'LRG_321p10'}, {0: 'LRG_321', 1: '+', 3: 'LRG_321t3', 4: 'LRG_321p3'}, {0: 'LRG_321', 1: '+', 3: 'LRG_321t2', 4: 'LRG_321p2'}, {0: 'LRG_321', 1: '+', 3: 'LRG_321t1', 4: 'LRG_321p1'}]")

    def test_UniProtAdapter(self):
        self.assertWarnings(DeprecationWarning, UniProtDB)

    def test_RefSeqAdapter(self):
        self.assertWarnings(DeprecationWarning, RefSeqAdapter,"1","2","3","4")