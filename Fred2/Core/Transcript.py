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
    """A Transcript is the mRNA sequence containing no or several :class:`Fred2.Core.Variant.Variant`.

    .. note::

        For accessing and manipulating the sequence see also :mod:`Bio.Seq.Seq`
        (from Biopython)

    :param str gene_id: Genome ID
    :param str transcript_id: :class:`~Fred2.Core.Transcript.Transcript` RefSeqID
    :param vars: Dict of :class:`Fred2.Core.Variant.Variant` for specific positions in the
                 :class:`~Fred2.Core.Transcript.Transcript`. key=position, value=Variant
    :type vars: dict(int,:class:`Fred2.Core.Variant.Variant`)
    """
    newid = itertools.count().__next__ #this is evil

    def __init__(self, seq, gene_id="unknown", transcript_id=None, vars=None):
        """
        :param str gene_id: Genome ID
        :param str transcript_id: :class:`~Fred2.Core.Transcript.Transcript` RefSeqID
        :param str seq: :class:`~Fred2.Core.Transcript.Transcript` sequence
        :param vars: A dict of :class:`~Fred2.Core.Transcript.Transcript` position to :class:`Fred2.Core.Variant.Variant`
                     that is specific to the :class:`~Fred2.Core.Transcript.Transcript`
        :type vars: dict(int,:class:`Fred2.Core.Variant.Variant`)
        """
        MetadataLogger.__init__(self)
        Seq.__init__(self, seq.upper(), generic_rna)
        self.gene_id = gene_id
        self.transcript_id = Transcript.newid() if transcript_id is None else transcript_id
        self.vars = dict() if vars is None else vars

    def __getitem__(self, index):
        """
        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)

         Allows only simple slicing (i.e. start < stop)

        :param int index: position 
        :returns: A :class:`~Fred2.Core.Transcript.Transcript` consisting of the single letter at position :attr:`index`
                  or a new sliced :class:`~Fred2.Core.Transcript.Transcript` (following :mod:`Bio.Seq.Seq` definition)
        :rtype: :class:`~Fred2.Core.Transcript.Transcript`
        """
        if isinstance(index, int):
            #Return a single letter as a string
            return str(self)[index]
        else:
            start, stop, step = index.indices(len(self))
            if start > stop:
                raise ValueError("start has to be greater than stop")

            if index.step:
                slice = set(range(start, step, stop))
            else:
                slice = set(range(start, stop))

            _vars = {}
            _fs = {}
            shift = 0
            #collect also all frame shift variants that are not canceled out
            for pos, v in sorted(self.vars.items()):
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
            t = Transcript(seq, gene_id=self.gene_id, transcript_id=trans_id, vars=_vars)
            return t

    def __repr__(self):
        lines = ["TRANSCRIPT: %s" % self.transcript_id]
        # get all variants:
        lines += ["VARIANTS:"]
        for vpos, var in self.vars.items():
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