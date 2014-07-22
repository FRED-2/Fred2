# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'walzer', 'haegele', 'schubert', 'szolek'

import string
import csv
import re
import urllib2
import warnings
import logging
import operator
from itertools import izip
from operator import itemgetter, attrgetter

import vcf

from Core.Variant import Variant, MutationSyntax
from Core.Transcript import Transcript, TranscriptSet
from Core.Allele import Allele, AlleleSet
from Core.Peptide import Peptide, PeptideSet
from RefSeqDB import RefSeqDB


def check_min_req_GSvar(row):
    """
    checking the presence of mandatory columns
    :param row: dictionary of a GSvar row
    :return: boolean, True if min req met
    """
    if ("#chr" in row.keys()
                and "start" in row.keys()
                and "end" in row.keys()
                and "ref" in row.keys()
                and "obs" in row.keys()
                and "coding_and_splicing_details" in row.keys()
                and "variant_details" in row.keys()):
        return True
    return False


def read_GSvar(filename, sample_id, just_dict=False):
    """
    reads GSvar and tsv files (tab sep files in context of genetic variants), omitting and warning about rows missing
    mandatory columns
    :param filename: /path/to/file
    :param sample_id:
    :param just_dict:
    :return: list of dictionaries representing valid variations
    """
    # TODO param type check
    ld_var = list()
    with open(filename, 'rb') as tsvfile:
        tsvreader = csv.DictReader((row for row in tsvfile if not row.startswith('##')), dialect='excel-tab')
        for row in tsvreader:
            if not check_min_req_GSvar(row):
                logging.warn("read_GSvar: Omitted row! Mandatory columns not present in: \n"+str(row))
                #https://docs.python.org/2/library/warnings.html
            else:
                ld_var.append(row)

    #test = sorted(ld_var,key=itemgetter('#chr','start'))
    if just_dict:
        return ld_var

    var_list = list()
    for v in ld_var:
        var_list.append(Variant(v["#chr"], v["start"], v["end"], v["ref"], v["obs"], sample_id, v))
    return var_list


def parse_annotation(annotations, details=None):
    """
    parses an annotation string (http://www.hgvs.org/mutnomen/quickref.html)
        - e.g. from ANNOVAR: (string) ANK3:NM_001204403:exon10:c.939delG:p.M313fs,ANK3:NM_020987:exon9:c.957delG:p.M319fs,ANK3:NM_001204404:exon9:c.906delG:p.M302fs,
        - e.g. cosmic: (list of dict) ITGAV	ENST00000261023	71212	c.1502C>A	p.S501Y	...
    :param annotations: the ANNOVAR annotation string
    :param details: the ANNOVAR variant details
    :return:
    """
    mut_syn = dict()
    #gss - old ANNOVAR stuff
    # HAT1:NM_003642:exon11:c.A1208T:p.Q403L,
    # NBPF10:NM_001039703:exon4:p.A179D:GCT->GAT:+,

    #gsvar - ANNOVAR uses UCSC?
    # RANBP2:NM_006267:exon14:c.A1954G:p.I652V,
    # ANAPC1:NM_022662:exon12:c.C1393T:p.Q465X,
    # NCBP1:NM_002486:exon9:c.913_914insG:p.G305fs,
    # ANK3:NM_001204403:exon10:c.939delG:p.M313fs,ANK3:NM_020987:exon9:c.957delG:p.M319fs,ANK3:NM_001204404:exon9:c.906delG:p.M302fs,

    #cosmic - http://www.hgvs.org/mutnomen/ (http://www.hgvs.org/mutnomen/quickref.html)
    # Gene Name	Accession Number    CDS Mutation Syntax	AA Mutation Syntax
    # TP53	ENST00000269305	c.?	p.E286K
    # MAP3K6	NM_004672	c.2483C>T	p.T828I	...
    # ITGAV	ENST00000261023	c.1502C>A	p.S501Y	...
    # TP53	ENST00000269305	c.847_848insT	p.R283fs*23

    #maybe use HGVS 0.3.3dev-d13418a64c51 ??
    if isinstance(annotations, str):
        annotations = [x for x in annotations.split(',') if x]  # removes trailing list entry ''
        print annotations
        for a in annotations:
                    t, c, p = [None]*3
                    for i in a.split(':'):
                        if i.startswith(('NM_', 'ENST')):
                                t = i
                        elif i.startswith('c.'):
                                c = i
                        elif i.startswith('p.'):
                                p = i
                    if t and c and p:
                        mut_syn[t] = MutationSyntax(details, c, p)
                    else:
                        logging.warn('Insufficient data')

    elif isinstance(annotations, list) \
            and all(isinstance(x, dict)
                    and 'Accession Number' in x and 'CDS Mutation Syntax' in x and 'AA Mutation Syntax' in x
                    for x in annotations):
        for a in annotations:
            mut_syn[a['Accession Number']] = MutationSyntax(details, a['CDS Mutation Syntax'], a['AA Mutation Syntax'])
    else:
        logging.warn('Unusable input')
    print mut_syn
    return mut_syn


def read_vcf(filename):
    """
    reads vcf files and handles the pass through annovar producing list of dictionaries like readGSvar
    read in a vcf_file, save all available metadata in a dictionary (keys: metadata, infos, filters, formats, samples)
    as well as create and return a list of dicts including the data lines of the vcf file (each dict corresponds to a
    record using the column headers as keys)
    :param filename: /path/to/file
    :return: list of dictionaries
    """
    vcf_reader = vcf.Reader(open(filename, 'r'))
    metadata_dict = dict()
    metadata = vcf_reader.metadata
    infos = vcf_reader.infos
    filters = vcf_reader.filters
    formats = vcf_reader.formats
    samples = vcf_reader.samples
    metadata_dict["metadata"] = metadata
    metadata_dict["infos"] = infos
    metadata_dict["filters"] = filters
    metadata_dict["formats"] = formats
    metadata_dict["samples"] = samples

    list_records = []
    for record in vcf_reader:
        record_dict = dict()
        record_dict["CHROM"] = record.CHROM
        record_dict["POS"] = record.POS
        record_dict["ID"] = record.ID
        record_dict["REF"] = record.REF
        record_dict["ALT"] = record.ALT
        record_dict["QUAL"] = record.QUAL
        record_dict["FILTER"] = record.FILTER
        record_dict["INFO"] = record.INFO
        for i in range(len(vcf_reader.samples)):
            record_dict[str(vcf_reader.samples[i])] = record.samples[i]
        #output_line = '%i\t%i\t%i\t%s\t%s\t%s\n' % \
        #          (list_records)
        #print record.samples[0]
        #print record_dict
        list_records.append(record_dict)
    return list_records, metadata_dict


def import_allele_list(file):
    """
    reads a csv file containing alleles their population probabilities and binding affinity thresholds
    The content should look like this:

    MHC-name1,pop_prob,binding_thresh
    MHC-name2,pop_prob,binding_thresh

    :param file (String): the location of the allele list
    :return: alleleSet: an alleleSet containing those informations
    """

    alleleSet = AlleleSet()
    with open(file, "r") as f:

        for l in f:
            allele,prob,thr = l.strip().split(",")
            a = Allele(allele)
            a.log_metadata("prob",float(prob))
            a.log_metadata("bindingThresh",float(thr))
            alleleSet[allele]=a

    return alleleSet


def import_epitope_list(file):
    """
    reads a csv file containing epitopes, their binding affinity prediction for all alleles in the allele list
    The content should look like this:

    epitopes,protein,variation,prediction_tool,MHC-allele_1,MHC-allele_2,.... ,MHC-allele_n
    epitope_seq1,prot_name,var_name,pred_name,affinit_1,affinit_2,.......,affinit_n
    ...

    :param file (String): the location of the epitope list
    :return: peptideSet: a peptideSet containing all listed epitopes containing binding-affinity and variation information
    """
    peptideSet = PeptideSet()

    with open(file, "r") as f:

        alleles = f.readline().strip().split(",")[4:]

        for l in f:
            splits = l.strip().split(",")
            affinities = splits[4:]

            if splits[0] not in peptideSet:
                p = Peptide(splits[0], None)
            else:
                p = peptideSet[splits[0]][0]
            p.log_metadata("variation", splits[1]+"-"+splits[2])

            for allele, affinity in izip(alleles, affinities):
                p.log_metadata("affinity", (splits[3],allele,float(affinity)))

            peptideSet.add_peptide(p)

    return peptideSet


def write_annovar_file(list_records, filename):
    """
    write out annovar file using the list of dicts generated from reading in the vcf file
    :param list_records: list of dicts created in the function before
    :param filename: name and directory of file that is created in scope of this function
    """
    output_file = open(filename, 'w')
    for rec in list_records:
        end_pos = rec["POS"] + len(rec["REF"]) - 1
        obs = reduce(operator.add, rec["ALT"])

        output_line = '%s\t%i\t%i\t%s\t%s\t%s\n' % \
                  (rec["CHROM"], rec["POS"], end_pos, rec["REF"], obs , rec["INFO"])
        output_file.write(output_line)
    output_file.close()


def write_peptide_file(pepset, pepfile):
    with open(pepfile, 'w') as f:
        for pepseq in pepset:
            f.write(pepseq + '\n')


def create_transcripts(varset, combinations=None):
    """
    assumes the variations attributes are properly set
    :param varset:
    :param combinations:
    """
    refseq_db = RefSeqDB()
    assert isinstance(varset, list) \
        and all(isinstance(var, Variant) for var in varset), 'wrong input - should be list of Variant_s'

    genes = set([x.find_gene(refseq_db) for x in varset])

    soon_transcripts = list()
    for g in genes:
        soon_transcripts.append([var for var in varset if var.gene == g])

    for vl in soon_transcripts:
        #sort should be stable!
        vl.sort(key=attrgetter('chromosome'))
        vl.sort(key=attrgetter('start'))
        vl.sort(key=attrgetter('stop'))
        vl.sort(key=attrgetter('gene'))

    transcripts = list()
    for st in soon_transcripts:
        ids = refseq_db.get_variant_ids(st[0].gene)
        for id in ids:
            pi, ps, ti, ts = [None]*4
            if id['RefSeq mRNA [e.g. NM_001195597]']:
                ti = id['RefSeq mRNA [e.g. NM_001195597]']
                pi = id['RefSeq Protein ID [e.g. NP_001005353]']
            else:
                ti = id['RefSeq mRNA [e.g. NM_001195597]']
                pi = id['RefSeq Protein ID [e.g. NP_001005353]']
            ts = refseq_db.get_transcript_sequence(ti)
            ps = refseq_db.get_product_sequence(pi)
            tr = Transcript(ti, ts, pi, ps)
            for var in st:
                tr.add_variant(var)
            transcripts.append(tr)

    # t = """MERGKMAEAESLETAAEHERILREIESTDTACIGPTLRSVYDGEEHGRFMEKLETRIRNHDREIEKMCNFHYQGFVDSITELLKVRGEAQKLKNQVTDTNRKLQHEGKELVIAMEELKQCRLQQRNISATVDKLMLCLPVLEMYSKLRDQMKTKRHYPALKTLEHLEHTYLPQVSHYRFCKVMVDNIPKLREEIKDVSMSDLKDFLESIRKHSDKIGETAMKQAQQQRNLDNIVLQQPRIGSKRKSKKDAYIIFDTEIESTSPKSEQDSGILDVEDEEDDEEVPGAQDLVDFSPVYRCLHIYSVLGARETFENYYRKQRRKQARLVLQPPSNMHETLDGYRKYFNQIVGFFVVEDHILHTTQGLVNRAYIDELWEMALSKTIAALRTHSSYCSDPNLVLDLKNLIVLFADTLQVYGFPVNQLFDMLLEIRDQYSETLLKKWAGIFRNILDSDNYSPIPVTSEEMYKKVVGQFPFQDIELEKQPFPKKFPFSEFVPKVYNQIKEFIYACLKFSEDLHLSSTEVDDMIRKSTNLLLTRTLSNSLQNVIKRKNIGLTELVQIIINTTHLEKSCKYLEEFITNITNVLPETVHTTKLYGTTTFKDARHAAEEEIYTNLNQKIDQFLQLADYDWMTGDLGNKASDYLVDLIAFLRSTFAVFTHLPGKVAQTACMSACKHLATSLMQLLLEAEVRQLTLGALQQFNLDVRECEQFARSGPVPGFQEDTLQLAFIDLRQLLDLFIQWDWSTYLADYGQPNCKYLRVNPVTALTLLEKMKDTSRKNNMFAQFRKNERDKQKLIDTVAKQLRGLISSHHS"""
    # i = """NM_015189"""
    # vc = "c.A1842C"
    # pc = "p.L614F"
    #
    # from IO import IO
    # from Core import Variant, Transcript
    # V = Variant.Variant(1,234,234,'A','C',"whatsoever",{})
    # V.log_metadata('coding', 'EXOC6B:NM_015189:exon18:c.A1842C:p.L614F,')
    # V.find_coding()
    # T = Transcript.Transcript(i, 'tseq', 'pid' , t)
    # T.add_variant(V)
    # T.translate()



