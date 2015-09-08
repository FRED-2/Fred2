# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.Transcript
   :synopsis: Contains the Transcript Class.
.. moduleauthor:: brachvogel, schubert, walzer

"""

import itertools

from Bio.Seq import Seq
from Bio.Alphabet import generic_rna

from Fred2.Core.Protein import Protein
from Fred2.Core.Base import MetadataLogger


class Transcript(MetadataLogger, Seq):
    """A Transcript is the mRNA sequence containing at no or several variations.

    .. note::

        For accessing and manipulating the sequence see also :mod:`Bio.Seq.Seq`
        (from Biopython)

    :param str gene_id: Genome ID
    :param str transcript_id: Transcript RefSeqID
    :param dict(int,Variant) vars: Dict of Variants for specific positions in
                                   the transcript. key=position, value=Variant
    """
    newid = itertools.count().next #this is evil

    def __init__(self, _seq, _gene_id="unknown", _transcript_id=None, _vars=None):
        """
        :param str _gene_id: input genome ID
        :param str _transcript_id: input transcript RefSeqID
        :param str _seq: Transcript RefSeq sequence
        :param list(Variant) _vars: a list of Variants which internal positions are specific to the transcript.
        """
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq, generic_rna)
        self.gene_id = _gene_id
        self.transcript_id = Transcript.newid() if _transcript_id is None else _transcript_id
        #TODO: this is not what the doc string says:
        if _vars is not None:
            self.vars = {v.get_transcript_position(_transcript_id): v for v in _vars}

        else:
            self.vars = dict()

    def __getitem__(self, index):
        """
        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)

        :param int index: position 
        :returns: (Transcript) -- A Transcript consisting of the single
        letter at position :attr:`index`.
        """
        item = str(self)[index]
        return Transcript(self.vars, self.transcript_id, item)

    def __repr__(self):
        lines = ["TRANSCRIPT: %s" % self.transcript_id]
        # get all variants:
        lines += ["VARIANTS:"]
        for vpos, var in self.vars.iteritems():
            lines.append('\tpos %i: %s'%(vpos, var))
        lines += ["SEQUENCE: %s (mRNA)"%str(self)]
        return '\n\t'.join(lines)

    def __eq__(self, other):
        return str(self) == str(other)

    def __lt__(self, other):
        return str(self) <= str(other)

    def __ge__(self, other):
        return str(self) >= str(other)

    def __cmp__(self, other):
        return cmp(str(self), str(other))

    def __hash__(self):
        return hash(str(self))