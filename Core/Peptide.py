"""
.. module:: Peptide
   :synopsis: Contains the Peptide Class.
.. moduleauthor:: brachvogel

"""
__author__ = 'brachvogel'

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Fred2.Core.Base import MetadataLogger
#from collections import OrderedDict


class Peptide(MetadataLogger, Seq):

    def __init__(self, _seq): #, _proteins=None, _vars=None, _transcripts=None):
        """
        :param _seq: String sequence in one letter amino acid code
        :type _seq: str.
        :param _proteins: dict of transcript_Ids to protein instances.
        :type _proteins: {str: Fred2.Core.Protein.Protein}.
        :param _vars: dictionary for variants according to position.
        :type _vars: {int:Fred2.Core.Variant.Variant}.
        :param _transcripts: dictionary for transcripts according to their id.
        :type _transcripts: {str: Fred2.Core.Transcript.Transcript}.
        """
        
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq, IUPAC.IUPACProtein)
        self.proteins = {}
        self.vars = {}
        self.transcripts = {}


    def __getitem__(self, index):
        """

        :param index: (int) position in the peptide sequence
        :returns: Peptide -- A Peptide consisting of the single letter at
            position :attr:`index`.
        """
        item = self[index]
        new_pept = Peptide(item)
        new_pept.proteins = self.proteins
        new_pept.vars = self.vars
        new_pept.transcripts = self.transcripts
        return new_pept


    def __repr__(self):
        return "" # TODO "%s found on %s (%s)" % (str(self.seq), self.origin.origin.id, self.origin.origin.gene)


# class PeptideSet(MetadataLogger, OrderedDict):

#     def __init__(self, peptides=None):
#         MetadataLogger.__init__(self)
#         OrderedDict.__init__(self)
#         if peptides is not None:
#             for peptide in peptides:
#                 self.add_peptide(peptide)

#     def add_peptide(self, peptide):
#         pepseq = str(peptide.sequence)
#         if pepseq in self:
#             self[pepseq].append(peptide)
#         else:
#             self[pepseq] = [peptide]

#     def filter_sample_ids(self, sample_ids):
#         # only keeps peptides who contain at least one variant from given sample_ids.
#         for pepseq, peplist in self.items():
#             peplist = [p for p in peplist if any(
#                 (v.sample_id in sample_ids for v in p.get_metadata('variants')[0]))]
#             if peplist:
#                 self[pepseq] = peplist
#             else:
#                 del self[pepseq]

#     def merge(self, other_pepset):
#         for pepseq, peplist in other_pepset.iteritems():
#             if pepseq in self:
#                 self[pepseq].extend(peplist)
#             else:
#                 self[pepseq] = peplist
