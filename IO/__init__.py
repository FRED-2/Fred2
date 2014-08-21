# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'walzer', 'haegele', 'schubert', 'szolek'

import csv
import re
import logging
import operator
from itertools import izip
from operator import itemgetter, attrgetter

import vcf
from Bio import SeqIO

import Core
from Core.Variant import Variant, MutationSyntax
from Core.Base import AASequence, Score
from Core.Allele import Allele
from Core.Peptide import Peptide, PeptideSet


def write_peptide_file(pepset, file, length=None):
    assert all(isinstance(a, AASequence) for a in pepset), "No list of AASequence"
    with open(file, 'w') as f:
        if length:
            f.writelines([str(x.seq) + '\n' for x in Core.lengthrestrict_list(pepset, length)])
        else:
            f.writelines([str(x.seq) + '\n' for x in pepset])


def write_peptide_fasta(pepset, fasta_file, length=None):
    assert all(isinstance(a, AASequence) for a in pepset), "No list of AASequence"
    with open(fasta_file, 'w') as f:
        if length:
            SeqIO.write(Core.lengthrestrict_list(pepset, length), f, "fasta")
        else:
            SeqIO.write(pepset, f, "fasta")


def write_protein_fasta(proset, fasta_file):
    assert all(isinstance(a, AASequence) for a in proset), "No list of AASequence"
    with open(fasta_file, 'w') as f:
        SeqIO.write(proset, f, "fasta")


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
                and ("coding_and_splicing_details" in row.keys() or "coding" in row.keys())
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
                logging.warning("read_GSvar: Omitted row! Mandatory columns not present in: \n"+str(row))
                #https://docs.python.org/2/library/warnings.html
            else:
                ld_var.append(row)

    #test = sorted(ld_var,key=itemgetter('#chr','start'))
    if just_dict:
        return ld_var

    var_list = list()
    for v in ld_var:
        var_list.append(Variant(v["#chr"].lstrip('chr'), v["start"], v["end"], v["ref"], v["obs"], sample_id, v))
    return var_list


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
            allele = l.strip()
            #allele, prob, thr = l.strip().split(",")
            a = Allele(allele)
            #a.log_metadata("prob", float(prob))
            #a.log_metadata("bindingThresh", float(thr))
            alleleSet[allele] = a

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

