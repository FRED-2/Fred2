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

from Fred2.Core.Base import MetadataLogger
from Fred2.Core.Variant import VariationType


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
        :param dict(int,Variant) _vars: a dict of transcript position to Variant that is specific to the transcript.
        """
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq.upper(), generic_rna)
        self.gene_id = _gene_id
        self.transcript_id = Transcript.newid() if _transcript_id is None else _transcript_id
        #TODO: this is not what the doc string says:
        self.vars = dict() if _vars is None else _vars

    def __getitem__(self, index):
        """
        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)

         Allows only simple slicing (i.e. start < stop)

        :param int index: position 
        :returns: (Transcript) -- A Transcript consisting of the single
        letter at position :attr:`index` or a new sliced Transcript (following Bio.Seqs definition)
        """
        if isinstance(index, int):
            #Return a single letter as a string
            return str(self)[index]
        else:
            start, stop, step = index.indices(len(self))
            if start > stop:
                raise ValueError("start has to be greater than stop")

            if index.step:
                slice = set(xrange(start, step, stop))
            else:
                slice = set(xrange(start, stop))

            _vars = {}
            _fs = {}
            shift = 0
            #collect also all frame shift variants that are not canceled out
            for pos, v in sorted(self.vars.iteritems()):
                if pos < start:
                    if v.type in [VariationType.FSINS, VariationType.FSDEL]:
                        shift = (v.get_shift()+shift) % 3
                        if shift:
                            _fs.setdefault[pos-start] = v
                        else:
                            _fs.clear()
                if pos in slice:
                    _vars[pos-start] = v
            _vars.update(_fs)
            trans_id = self.transcript_id+":"+str(Transcript.newid())
            seq = str(self)[index]
            t = Transcript(seq, _gene_id=self.gene_id, _transcript_id=trans_id, _vars=_vars)
            return t

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