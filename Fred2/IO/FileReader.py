# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Reader
   :synopsis: Module handles reading of files. line reading, FASTA reading, annovar reading
.. moduleauthor:: brachvogel, schubert

"""

import os
import re
import warnings

import vcf
from Bio.SeqIO.FastaIO import SimpleFastaParser

from Fred2.Core.Peptide import Peptide
from Fred2.Core.Variant import Variant, VariationType, MutationSyntax


####################################
#       F A S T A  -  R E A D E R
####################################
def read_fasta(files, in_type=Peptide, id_position=1):
    """
    Generator function:

    Read a (couple of) peptide, protein or rna sequence from a FASTA file.
    User needs to specify the correct type of the underlying sequences. It can
    either be: Peptide, Protein or Transcript (for RNA).

    :param files: A (list) of file names to read in
    :in_type files: list(str) or str
    :param in_type: The type to read in
    :type in_type: :class:`~Fred2.Core.Peptide.Peptide` or :class:`~Fred2.Core.Transcript.Transcript`
                or :class:`~Fred2.Core.Protein.Protein`
    :param int id_position: the position of the id specified counted by |
    :returns: a list of the specified sequence type derived from the FASTA file sequences.
    :rtype: (list(:attr:`in_type`))
    :raises ValueError: if a file is not readable
    """

    if isinstance(files, str):
            files = [files]
    else:
            if any(not os.path.exists(f) for f in files):
                raise ValueError("Specified Files do not exist")

    collect = set()
    # open all specified files:
    for name in files:
        with open(name, 'r') as handle:
            # iterate over all FASTA entries:
            for _id, seq in SimpleFastaParser(handle):
                # generate element:
                try:
                    _id = _id.split("|")[id_position]
                except IndexError:
                   _id = _id

                try:
                    collect.add(in_type(seq.strip().upper(), transcript_id=_id))
                except TypeError:
                    collect.add(in_type(seq.strip().upper()))
    return list(collect)


####################################
#       L I N E  -  R E A D E R
####################################
def read_lines(files, in_type=Peptide):
    """
    Generator function:

    Read a sequence directly from a line. User needs to manually specify the 
    correct type of the underlying data. It can either be:
    Peptide, Protein or Transcript, Allele.

    :param files: a list of strings of absolute file names that are to be read.
    :in_type files: list(str) or str
    :param in_type: Possible in_type are :class:`~Fred2.Core.Peptide.Peptide`, :class:`~Fred2.Core.Protein.Protein`,
                 :class:`~Fred2.Core.Transcript.Transcript`, and :class:`~Fred2.Core.Allele.Allele`.
    :type in_type: :class:`~Fred2.Core.Peptide.Peptide` or :class:`~Fred2.Core.Protein.Protein` or
                :class:`~Fred2.Core.Transcript.Transcript` or :class:`~Fred2.Core.Allele.Allele`
    :returns: A list of the specified objects
    :rtype: (list(:attr:`in_type`))
    :raises IOError: if a file is not readable
    """

    if isinstance(files, str):
            files = [files]
    else:
            if any(not os.path.exists(f) for f in files):
                raise IOError("Specified Files do not exist")

    collect = set()
    #alternative to using strings is like: cf = getattr(Fred2.Core, "Protein"/"Peptide"/"Allele"/...all in core)
    for name in files:
        with open(name, 'r') as handle:
            # iterate over all lines:
            for line in handle:
                # generate element:
                collect.add(in_type(line.strip().upper()))

    return list(collect)


#####################################
#       A N N O V A R  -  R E A D E R
#####################################
def read_annovar_exonic(annovar_file, gene_filter=None, experimentalDesig=None):
    """
    Reads an gene-based ANNOVAR output file and generates :class:`~Fred2.Core.Variant.Variant` objects containing
    all annotated :class:`~Fred2.Core.Transcript.Transcript` ids an outputs a list :class:`~Fred2.Core.Variant.Variant`.

    :param str annovar_file: The path ot the ANNOVAR file
    :param list(str) gene_filter: A list of gene names of interest (only variants associated with these genes
                                  are generated)
    :return: List of :class:`~Fred2.Core.Variant.Variants fully annotated
    :rtype: list(:class:`~Fred2.Core.Variant.Variant`)
    """

    vars = []
    gene_filter = gene_filter if gene_filter is not None else []

    #fgd3:nm_001083536:exon6:c.g823a:p.v275i,fgd3:nm_001286993:exon6:c.g823a:p.v275i,fgd3:nm_033086:exon6:c.g823a:p.v275i
    #RE = re.compile("\w+:(\w+):exon\d+:c.(\D*)(\d+)_*(\d*)(\D\w*):p.\w+:\D*->\D*:(\D).*?,")
    #RE = re.compile("\w+:(\w+):exon\d+:c.(\D*)(\d+)_*(\d*)(\D\w*):p.(\D*)(\d+)_*(\d*)(\D\w*):(\D).*?,")
    RE = re.compile("((\w+):(\w+):exon\d+:c.\D*(\d+)\D\w*:p.\D*(\d+)\D\w*)")
    type_mapper = {('synonymous', 'snv'): VariationType.SNP,
                   ('nonsynonymous', 'snv'): VariationType.SNP,
                   ('stoploss', 'snv'): VariationType.SNP,
                   ('stopgain', 'snv'): VariationType.SNP,
                   ('nonframeshift', 'deletion'): VariationType.DEL,
                   ('frameshift', 'deletion'): VariationType.FSDEL,
                   ('nonframeshift', 'insertion'): VariationType.INS,
                   ('frameshift', 'insertion'): VariationType.FSINS}
    with open(annovar_file, "r") as f:
        for line in f:
            mut_id, mut_type, line, chrom, genome_start, genome_stop, ref, alt, zygos = [x.strip().lower() for x in line.split("\t")[:9]]
            #print ref, alt

            #test if its a intersting snp

            gene = line.split(":")[0].strip().upper()

            if gene not in gene_filter and len(gene_filter):
                continue

            if gene == "UNKNOWN":
                warnings.warn("Skipping UNKWON gene")
                continue

           # print "Debug ", gene, type.split(),mut_id
            #print "Debug ", line, RE.findall(line), type, zygos
            coding = {}
            for nm_id_pos in RE.findall(line):
                mutation_string, geneID, nm_id, trans_pos, prot_start = nm_id_pos
                #print "Debug ",nm_id_pos

                nm_id = nm_id.upper()
                _,_, _, trans_coding, prot_coding = mutation_string.split(":")
                #internal transcript and protein position start at 0!
                coding[nm_id] = MutationSyntax(nm_id, int(trans_pos)-1, int(prot_start)-1, trans_coding, prot_coding,
                                               geneID=geneID.upper())

            ty = tuple(mut_type.split())

            vars.append(
                Variant(mut_id, type_mapper.get(ty, VariationType.UNKNOWN), chrom, int(genome_start), ref.upper(),
                        alt.upper(), coding, zygos == "hom", ty[0] == "synonymous",
                        experimentalDesign=experimentalDesig))
    return vars


#####################################
#       V C F  -  R E A D E R
#####################################
def read_vcf(vcf_file, gene_filter=None, experimentalDesig=None):
    """
    Reads an vcf v4.0 or 4.1 file and generates :class:`~Fred2.Core.Variant.Variant` objects containing
    all annotated :class:`~Fred2.Core.Transcript.Transcript` ids an outputs a list :class:`~Fred2.Core.Variant.Variant`.
    Only the following variants are considered by the reader where synonymous labeled variants will not be integrated into any variant:
    filter_variants = ['missense_variant', 'frameshift_variant', 'stop_gained', 'missense_variant&splice_region_variant', "synonymous_variant", "inframe_deletion", "inframe_insertion"]

    :param str vcf_file: The path ot the vcf file
    :param list(str) gene_filter: A list of gene names of interest (only variants associated with these genes
                                  are generated)
    :return: List of :class:`~Fred2.Core.Variant.Variants fully annotated
    :rtype: Tuple of (list(:class:`~Fred2.Core.Variant.Variant`), list(transcript_ids)
    """
    vl = list()
    with open(vcf_file, 'rb') as tsvfile:
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        vl = [r for r in vcf_reader]

    list_vars = []
    transcript_ids = []

    genotye_dict = {"het": False, "hom": True, "ref": True}

    for num, record in enumerate(vl):
        c = record.CHROM.strip('chr')  # chrom
        p = record.POS - 1  # vcf is 1-based & FRED2 0-based
        variation_dbid = record.ID  # e.g. rs0123
        r = str(record.REF)  # reference nuc (seq)
        v_list = record.ALT  # list of variants
        q = record.QUAL  # ?
        f = record.FILTER  # empty if PASS, content otherwise
        # I guess we shouldn't expect that keyword to be there ?!
        #z = record.INFO['SOMATIC'] #if true somatic

        vt = VariationType.UNKNOWN
        if record.is_snp:
            vt = VariationType.SNP
        elif record.is_indel:
            if len(v_list)%3 == 0:  # no frameshift
                if record.is_deletion:
                    vt = VariationType.DEL
                else:
                    vt = VariationType.INS
            else:  # frameshift
                if record.is_deletion:
                    vt = VariationType.FSDEL
                else:
                    vt = VariationType.FSINS
        gene = None

        # WHICH VARIANTS TO FILTER ?
        filter_variants = ['missense_variant', 'frameshift_variant', 'stop_gained', 'missense_variant&splice_region_variant', "synonymous_variant", "inframe_deletion", "inframe_insertion"]

        for alt in v_list:
            isHomozygous = False
            if 'HOM' in record.INFO:
                #TODO set by AF & FILTER as soon as available
                isHomozygous = record.INFO['HOM'] == 1
            elif 'SGT' in record.INFO:
                zygosity = record.INFO['SGT'].split("->")[1]
                if zygosity in genotye_dict:
                    isHomozygous = genotye_dict[zygosity]
                else:
                    if zygosity[0] == zygosity[1]:
                        isHomozygous = True
                    else:
                        isHomozygous = False
            else:
                for sample in record.samples:
                    if 'GT' in sample.data:
                        isHomozygous = sample.data['GT'] == '1/1'

            if "ANN" in record.INFO and record.INFO['ANN']:
                isSynonymous = False
                coding = dict()
                for annraw in record.INFO['ANN']:  # for each ANN only add a new coding! see GSvar
                    annots = annraw.split('|')

                    obs, a_mut_type, impact, a_gene, a_gene_id, feature_type, transcript_id, exon, tot_exon, trans_coding, prot_coding, cdna, cds, aa, distance, warns = annots

                    if a_mut_type in filter_variants:
                        tpos = 0
                        ppos = 0

                        # get cds/protein positions and convert mutation syntax to FRED2 format
                        if trans_coding != '':
                            positions = re.findall(r'\d+', trans_coding)
                            ppos = int(positions[0]) - 1

                        if prot_coding != '':
                            positions = re.findall(r'\d+', prot_coding)
                            tpos = int(positions[0]) - 1

                        isSynonymous = (a_mut_type == "synonymous_variant")

                        #rather take gene_id than gene name
                        gene = a_gene_id

                        #REFSEQ specific ? Do have to split because of biomart ?
                        transcript_id = transcript_id.split(".")[0]

                        #TODO vcf are not REFSEQ only

                        #coding string not parsed anyway ? just use the one given by SnpEff
                        coding[transcript_id] = MutationSyntax(transcript_id, ppos, tpos, trans_coding, prot_coding)
                        transcript_ids.append(transcript_id)

                if coding and not isSynonymous:
                    if vt == VariationType.SNP:
                        pos, reference, alternative = p, str(r), str(alt)
                    elif vt == VariationType.DEL or vt == VariationType.FSDEL:
                        if alt != '-':
                            pos, reference, alternative = p + len(alt), r[len(alt):], '-'
                        else:
                            pos, reference, alternative = p, str(r), str(alt)
                    elif vt == VariationType.INS or vt == VariationType.FSINS:
                        if r != '-':
                            if alt != '-':
                                pos, reference, alternative = p + len(r), '-', str(alt)[len(r):]
                            else:
                                pos, reference, alternative = p + len(r), '-', str(alt)
                        else:
                            pos, reference, alternative = p, str(r), str(alt)

                    var = Variant("line" + str(num), vt, c, pos, reference, alternative, coding, isHomozygous, isSynonymous, experimentalDesign=experimentalDesig)
                    var.gene = gene
                    var.log_metadata("vardbid", variation_dbid)
                    list_vars.append(var)

            else:
                warnings.warn("Skipping unannotated variant", UserWarning)

    return list_vars, transcript_ids