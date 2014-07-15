# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
import string
import csv
import re
import urllib2
import warnings
import logging
import operator
from itertools import izip

import vcf

from Core.Variant import Variant
from Core.Transcript import Transcript, TranscriptSet
from Core.Allele import Allele, AlleleSet
from Core.Peptide import Peptide, PeptideSet

__author__ = 'walzer', 'haegele', 'schubert', 'szolek'


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


def read_GSvar(filename):
    """
    reads GSvar and tsv files (tab sep files in context of genetic variants), omitting and warning about rows missing
    mandatory columns
    :param filename: /path/to/file
    :return: list of dictionaries representing valid variations
    """
    ld_var = list()
    with open(filename, 'rb') as tsvfile:
        tsvreader = csv.DictReader((row for row in tsvfile if not row.startswith('##')), dialect='excel-tab')
        for row in tsvreader:
            if not check_min_req_GSvar(row):
                logging.warn("readGSvar: Omitted row! Mandatory columns not present in: \n"+str(row))
                #https://docs.python.org/2/library/warnings.html
            else:
                ld_var.append(row)
    #test = sorted(ld_var,key=itemgetter('#chr','start'))
    return ld_var


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


def parse_annovar_annotation(annotation, details):
    """
    parses a annotation string from ANNOVAR, needs the variant details of annovar as well
    as the annotation string has different formatting for the variants.
    :param annotation: the ANNOVAR annotation string
    :param details: the ANNOVAR variant details
    :return:
    """
    return None


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


def import_annovar_tsv(annovar_file, sample_id, refseq_lookup_handle, max_lines=100000):

    complement = string.maketrans('atgcATGC', 'tacgTACG')
    ts = TranscriptSet()

    data = csv.reader(open(annovar_file), delimiter='\t')
    fields = data.next()
    for jj, onerow in enumerate(data):

        if jj > max_lines:
            break  # limit to X lines for quick testing

        row = dict(zip(fields, onerow))
        genomic_pos = (row['#chr'], row['start'], row['end'])

        # empty sequences show up as a single '-' in the ref and obs columns. Delete that dash.
        seq_obs = row['obs'].strip('-')
        seq_ref = row['ref'].strip('-')

        if row['variant'] == 'SNV':
            variant_type = 'SNV'
            assert len(seq_obs) == 1 and len(seq_ref) == 1, 'SNV observed and reference have improper lengths'
        elif row['variant'] == 'INDEL':
            if seq_obs:  # insertion
                assert seq_ref == '', 'insertion INDEL has a non-empty reference sequence'
                variant_type = 'FSI' if len(seq_obs) % 3 else 'INS'  # frameshift or not
            else:
                assert seq_obs == '', 'deletion INDEL has a non-empty observed sequence'
                variant_type = 'FSD' if len(seq_ref) % 3 else 'DEL'
        else:
            print 'variant type not SNV nor INDEL. What the hell?'

        # reverse complements for reverse strand transcripts
        seq_obs_rc = seq_obs[::-1].translate(complement)
        seq_ref_rc = seq_ref[::-1].translate(complement)

        x1 = Variant(genomic_pos, variant_type, seq_obs, 'N', sample_id)
        x1_rc = Variant(genomic_pos, variant_type, seq_obs_rc, 'N', sample_id)

        # if the variant is homozygous we know that the reference sequence cannot be present
        # in our sample so it's from HG19. Important to not mix it with our sample_ids.
        # Homozygous indels really bother me though. Why are there so many? What to do with them? TODO
        x2 = Variant(genomic_pos, 'REF', seq_ref, 'N', 'hg19'
            if row['genotype'] == 'hom' else (sample_id + '_ref'))  # TODO: I don't like it
        # But how should we represent the non-variant base in het. tumor mutations?
        x2_rc = Variant(genomic_pos, 'REF', seq_ref_rc, 'N', 'hg19'
            if row['genotype'] == 'hom' else (sample_id + '_ref'))  # TODO: I don't like it

        x1.log_metadata('genotype', row['genotype'])
        x2.log_metadata('genotype', row['genotype'])
        x1_rc.log_metadata('genotype', row['genotype'])
        x2_rc.log_metadata('genotype', row['genotype'])

        x1.log_metadata('original_row', (sample_id + '_' + str(jj), onerow))  # TODO: don't
        x2.log_metadata('original_row', (sample_id + '_' + str(jj), onerow))
        x1_rc.log_metadata('original_row', (sample_id + '_' + str(jj), onerow))
        x2_rc.log_metadata('original_row', (sample_id + '_' + str(jj), onerow))

        vs = VariantSet(genomic_pos, [x1, x2])
        vs_rc = VariantSet(genomic_pos, [x1_rc, x2_rc])

        transcripts = row['coding'].rstrip(',').split(',')
        for tr in transcripts:
            try:
                gene, tr_id, _, mutation_nuc, mutation_aa = tr.split(':')
                if tr_id in ts:
                    transcript = ts[tr_id]
                else:
                    try:
                        transcript = Transcript(tr_id, lookup_fn=refseq_lookup_handle)
                    except AssertionError as aerr:  # TODO: replace this to a proper exception.
                        print aerr
                        continue

                    ts.add_transcript(transcript)

                # a single number for SNV-s, two consecutive numbers for INS-es and two unconstrained
                # numbers for deletions. We turn them into a (start, stop) tuple, indexed from zero
                # and referring to positions on the reference sequence. In the following form:
                # SNVs: (x, x+1) - pinpointing the affected nucleotide
                # INSs: (x, x)   - 0 affected bases in the reference insertion spot is well defined
                # DELs: (x, y)   - deletion starting somewhere and ending somewhere
                # ...making it possible to define substituted slices in a concise way.
                positions = map(int, re.findall(r'[0-9_]+', mutation_nuc)[0].split('_'))
                positions.append(positions[-1])  # padding for single nucleotide deletion cases
                if variant_type == 'SNV':
                    vpos = (positions[0] - 1, positions[0])
                elif variant_type in ('INS', 'FSI'):
                    vpos = (positions[0], positions[0])
                elif variant_type in ('DEL', 'FSD'):
                    vpos = (positions[0] - 1, positions[1])

                # TODO: vpos is still relative to the CDS start position. It should probably be
                # relative to the transcript start position.
                transcript.add_variantset(vpos, vs if transcript.strand == '+' else vs_rc)

            except ValueError as ee:
                ee == ee  # don't need it now
                #print 'unparseable stuff: ', row['gene']
                continue
    return ts


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

