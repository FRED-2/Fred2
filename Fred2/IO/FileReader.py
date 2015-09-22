# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Reader
   :synopsis: Module handles reading of files. line reading, FASTA reading, annovar reading
.. moduleauthor:: brachvogel, schubert

"""

import warnings
import os
import re

from Bio.SeqIO.FastaIO import SimpleFastaParser

from Fred2.Core.Peptide import Peptide
from Fred2.Core.Variant import Variant, VariationType, MutationSyntax

####################################
#       F A S T A  -  R E A D E R
####################################
def read_fasta(files, type=Peptide, id_position=1):
    """
    Generator function:

    Read a (couple of) peptide, protein or rna sequence from a FASTA file.
    User needs to specify the correct type of the underlying sequences. It can
    either be: Peptide, Protein or Transcript (for RNA).

    :param List(str)/str files: A (List) of file names to read in
    :param Peptide/Transcript/Protein type: The type to read in
    :param int id_position: the position of the id specified counted by |
    :returns: (list(SequenceType)) -- a list of the specified sequence type
              derived from the FASTA file sequences.
    :exception ValueError: if a file is not readable
    """

    if isinstance(files, basestring):
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
                except KeyError:
                   _id = _id

                try:
                    collect.add(type(seq.strip().upper(), _transcript_id=_id))
                except TypeError:
                    collect.add(type(seq.strip().upper()))
    return list(collect)



####################################
#       L I N E  -  R E A D E R
####################################
def read_lines(files, type=Peptide):
    """
    Generator function:

    Read a sequence directly from a line. User needs to manually specify the 
    correct type of the underlying data. It can either be:
    Peptide, Protein or Transcript, Allele.

    :param list(str)/str files: a list of strings of absolute file names of the FASTA files
                  that are to be read. Give as: 
                  get_sequence(*["path/name1", "path/name2"]).
                  This field is required!
    :param Peptide/Protein/Transcript/Allele type: optional. Use get_sequence(*["path/name1", "path/name2",
                     type="Protein"). Possible types are Peptide, Protein, Transcript and Allele.
    :returns: (list(type)) -- a list of the specified objects
    :exception: if a file is not readable
    """

    if isinstance(files, basestring):
            files = [files]
    else:
            if any(not os.path.exists(f) for f in files):
                raise IOError("Specified Files do not exist")

    collect = set()
    for name in files:
        with open(name, 'r') as handle:
            # iterate over all lines:
            for line in handle:
                # generate element:
                collect.add(type(line.strip().upper()))

    return list(collect)


#####################################
#       A N N O V A R  -  R E A D E R
#####################################
def read_annovar_exonic(annovar_file, gene_filter=None, experimentalDesig=None):
    """
    Reads an gene-based ANNOVAR output file and generates Variant objects containing
    all annotated transcript ids an outputs a list variants.

    :param str annovar_file: The path ot the ANNOVAR file
    :param list(str) gene_filter: A list of gene names of interest (only variants associated
                                  with these genes are generated)
    :return: list(AVariant) -- List of variants fully annotated
    """

    vars = []
    gene_filter = gene_filter if gene_filter is not None else []

    #fgd3:nm_001083536:exon6:c.g823a:p.v275i,fgd3:nm_001286993:exon6:c.g823a:p.v275i,fgd3:nm_033086:exon6:c.g823a:p.v275i
    #RE = re.compile("\w+:(\w+):exon\d+:c.(\D*)(\d+)_*(\d*)(\D\w*):p.\w+:\D*->\D*:(\D).*?,")
    #RE = re.compile("\w+:(\w+):exon\d+:c.(\D*)(\d+)_*(\d*)(\D\w*):p.(\D*)(\d+)_*(\d*)(\D\w*):(\D).*?,")
    RE = re.compile("((\w+):exon\d+:c.\D*(\d+)\D\w*:p.\D*(\d+)\D\w*)")
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
            id, type, line, chrom, genome_start, genome_stop, ref, alt, zygos = map(lambda x: x.strip().lower(),
                                                                                    line.split("\t")[:9])
            #print ref, alt

            #test if its a intersting snp

            gene = line.split(":")[0].strip().upper()

            if gene not in gene_filter and len(gene_filter):
                continue

            if gene == "UNKNOWN":
                warnings.warn("Skipping UNKWON gene")
                continue

           # print "Debug ", gene, type.split(),id
            #print "Debug ", line, RE.findall(line), type, zygos
            coding = {}
            for nm_id_pos in RE.findall(line):
                mutation_string, nm_id, trans_pos, prot_start = nm_id_pos
                #print "Debug ",nm_id_pos

                nm_id = nm_id.upper()
                _, _, trans_coding, prot_coding = mutation_string.split(":")
                #internal transcript and protein position start at 0!
                coding[nm_id] = MutationSyntax(nm_id, int(trans_pos)-1, int(prot_start)-1, trans_coding, prot_coding)

            ty = tuple(type.split())

            vars.append(
                Variant(id, type_mapper.get(ty, VariationType.UNKNOWN), chrom, int(genome_start), ref.upper(),
                        alt.upper(), coding, zygos == "hom", ty[0] == "synonymous",
                        experimentalDesign=experimentalDesig))
    return vars
