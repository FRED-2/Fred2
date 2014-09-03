# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: FastA
   :synopsis: Module handles fasta file IO
.. moduleauthor:: brachvogel

"""

from Bio.SeqIO.FastaIO import SimpleFastaParser

from Fred2.Core.Peptide import Peptide
from Fred2.Core.Protein import Protein
from Fred2.Core.Transcript import Transcript

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
    type = kwargs.get("type", "Peptide")

    print type
    def _create_single_sequ(sequ, id):
        if type == "Peptide":
            return Peptide(sequ)
    
        elif type == "Protein":
            return Protein(sequ, "unknown", id)
    
        elif type == "Transcript":
            return Transcript("unkown", id, sequ)
    
        else:
            raise ValueError("An invalid sequence object type was specified for\
parsing a FASTA file. Type was %s but allowed types are: 'Peptide', 'Protein' \
and 'Transcript'."%type)


    ####################################
    res = [] # result list

    # open all specified files:
    for name in argv:
        with open(name, 'r') as handle:
            # iterate over all FASTA entries:
            for id, seq in SimpleFastaParser(handle):
                # generate element:
                seq = seq.upper()
                res.append(_create_single_sequ(seq, id))
    return res

# An example:
# if __name__ == '__main__':
#     root_dir = "/Users/pbrach/files/projects/Hiwi-Kohlbacher-2014/Galaxy/develop_fred2/"
#     files = [root_dir + "peptides_1to2.fasta",
#              root_dir + "peptides_3to4.fasta",
#              root_dir + "peptides_5to6.fasta",
#              root_dir + "peptides_7to8.fasta"]

#     single_file_lower = root_dir + "peptides_1to8_lower.fasta"
#     get_sequence(*files)
#     get_sequence(single_file_lower, type="Peptide")







