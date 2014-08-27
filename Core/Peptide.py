"""
.. module:: Peptide
   :synopsis: Contains the Peptide Class.
.. moduleauthor:: brachvogel

"""

from Base import MetadataLogger
from collections import OrderedDict


class Peptide(MetadataLogger):

    def __init__(self, sequence, origin_protein):
        MetadataLogger.__init__(self)
        self.sequence = sequence  # TODO: force Seq and protein_alphabet
        # self.length = len(sequence)
        self.origin = origin_protein

        # TODO: sample_id and variants found on peptide would be nice to know

        # external modifier class to handle scoring files and assigning
        # objects a score attribute

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        return "%s found on %s (%s)" % (str(self.sequence), self.origin.origin.id, self.origin.origin.gene)


class PeptideSet(MetadataLogger, OrderedDict):

    def __init__(self, peptides=None):
        MetadataLogger.__init__(self)
        OrderedDict.__init__(self)
        if peptides is not None:
            for peptide in peptides:
                self.add_peptide(peptide)

    def add_peptide(self, peptide):
        pepseq = str(peptide.sequence)
        if pepseq in self:
            self[pepseq].append(peptide)
        else:
            self[pepseq] = [peptide]

    def filter_sample_ids(self, sample_ids):
        # only keeps peptides who contain at least one variant from given sample_ids.
        for pepseq, peplist in self.items():
            peplist = [p for p in peplist if any(
                (v.sample_id in sample_ids for v in p.get_metadata('variants')[0]))]
            if peplist:
                self[pepseq] = peplist
            else:
                del self[pepseq]

    def merge(self, other_pepset):
        for pepseq, peplist in other_pepset.iteritems():
            if pepseq in self:
                self[pepseq].extend(peplist)
            else:
                self[pepseq] = peplist
