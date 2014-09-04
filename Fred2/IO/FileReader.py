# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Reader
   :synopsis: Module handles reading of files. line reading and FASTA reading
.. moduleauthor:: brachvogel

"""

from Bio.SeqIO.FastaIO import SimpleFastaParser

from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Protein import Protein
from Fred2.Core.Transcript import Transcript

def __check_type(type, allowed=['Peptide', 'Protein','Transcript']):
    """
    :param str type: the wrong type
    """
    if not type in allowed:
        raise ValueError("An invalid sequence object type was specified for \
parsing a FASTA file. Type was %s but allowed types are: %s."%(type, allowed))



def __create_single_sequ(sequ, id, type):
    if type == "Peptide":
        return Peptide(sequ)

    elif type == "Protein":
        return Protein(sequ, "unknown", id)

    elif type == "Transcript":
        return Transcript("unkown", id, sequ)

    elif type == "Allele":
        return Allele(sequ)

####################################
#       F A S T A  -  R E A D E R
####################################
def get_sequence(*argv, **kwargs):
    """
    Read a (couple of) peptide, protein or rna sequence from a FASTA file.
    User needs to specify the correct type of the underlying sequences. It can
    either be: Peptide, Protein or Transcript (for RNA).

    :param *argv: a string list of absolute file names of the FASTA files to be
                  read. Give as: get_sequence(*["path/name1", "path/name2"]).
                  Alternatively 
                  This field is required!
    :param **kwargs: optional. Use get_sequence(*["path/name1", "path/name2",
                     type="Protein"). Possible types are Peptide', 'Protein'
                     and 'Transcript'.
    :returns: (list(SequenceType)) -- a list of the specified sequence type
              derived from the FASTA file sequences.
    """
    _type = kwargs.get("type", "Peptide")

    __check_type(_type)

    collect = {}
    # open all specified files:
    for name in argv:
        with open(name, 'r') as handle:
            # iterate over all FASTA entries:
            for _id, seq in SimpleFastaParser(handle):
                # generate element:
                collect[seq.upper()] = _id
    
    return [__create_single_sequ(seq, _id, _type) \
            for seq, _id in collect.items()]

####################################
#       L I N E  -  R E A D E R
####################################
def read_lines(*argv, **kwargs):
    """
    Read a sequence directly from a line.
    User needs to specify the correct type of the underlying sequences. It can
    either be: Peptide, Protein or Transcript (for RNA).

    :param *argv: a string list of absolute file names of the FASTA files to be
                  read. Give as: get_sequence(*["path/name1", "path/name2"]).
                  Alternatively 
                  This field is required!
    :param **kwargs: optional. Use get_sequence(*["path/name1", "path/name2",
                     type="Protein"). Possible types are Peptide', 'Protein'
                     and 'Transcript'.
    :returns: (list(SequenceType)) -- a list of the specified sequence type
              derived from the FASTA file sequences.
    """
    _type = kwargs.get("type", "Peptide")

    __check_type(_type, ['Peptide', 'Protein', 'Transcript', 'Allele'])

    collect = set()
    for name in argv:
        with open(name, 'r') as handle:
            # iterate over all lines:
            for line in handle:
                # generate element:
                collect.add(line.strip().upper())

    return [__create_single_sequ(seq, "generic_id_"+str(_id), _type) \
            for _id, seq in enumerate(collect)]




# if __name__ == '__main__':
#     root_dir = "/Users/pbrach/files/projects/Hiwi-Kohlbacher-2014/Galaxy/develop_fred2/"
#     files = [root_dir + "peptides_1to2.fasta",
#              root_dir + "peptides_3to4.fasta",
#              root_dir + "peptides_5to6.fasta",
#              root_dir + "peptides_7to8.fasta"]

#     single_file_lower = root_dir + "peptides_1to8_lower.fasta"

#     objects = read_lines(single_file_lower, type="Transcript")
#     for elem in objects:
#         print elem.__repr__()
