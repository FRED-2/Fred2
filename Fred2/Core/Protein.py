# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.Protein
   :synopsis: Contains the Protein class

   :Note: All internal indices start at 0!

.. moduleauthor:: schubert, brachvogel, walzer

"""
import itertools

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Fred2.Core.Base import MetadataLogger
from Fred2.Core.Variant import VariationType


class Protein(MetadataLogger, Seq):
    """
    :class:`~Fred2.Core.Protein.Protein` corresponding to exactly one transcript.

    .. note::

        For accessing and manipulating the sequence see also :mod:`Bio.Seq.Seq`
        (from Biopython)

    """
    newid = itertools.count().__next__ #this is evil and has no other purpose? it does not help that there may be more than one protein from one transcript - due to variants

    def __init__(self, _seq, gene_id="unknown", transcript_id=None, orig_transcript=None, vars=None):
        """
        :param str _seq: String of an IUPACProtein alphabet, representing the protein
        :param str gene_id: ID of the genome the protein originated from
        :param str transcript_id: ID of the transcript the protein originated from
        :param orig_transcript: Reference to the originating transcript object
        :type orig_transcript: :class:`~Fred2.Core.Transcript.Transcript`
        :param vars: Nonsynonymous variants that are associated with the protein. key=position within protein,
                     value=list of variants at that pos
        :type vars: dict(int,list(:class:`~Fred2.Core.Variant.Variant`))
        """
        # Init parent type:
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq.upper(), IUPAC.IUPACProtein)
        # Init own member:
        if vars is None:
            self.vars = dict()
        else:
            self.vars = vars  # {prot-position: list(variant)}
        self.orig_transcript = orig_transcript
        self.transcript_id = "Protein_%i"%Protein.newid() if transcript_id is None else transcript_id
        self.gene_id = gene_id

    def __getitem__(self, index):
        """

        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)

        Returns a single letter or a sliced :class:`~Fred2.Core.Protein.Protein` (when given a slice).
        The sliced :class:`~Fred2.Core.Protein.Protein` does not reference to a
        :class:`~Fred2.Core.Transcript.Transcript` object

        Allows only simple slicing (i.e. start < stop)

        :param int/Slice index: position within the primary sequence or a slice
        :returns: A single letter at position :attr:`index` or a sliced :class:`~Fred2.Core.Protein.Protein` with
                  adjusted variant positions
        :rtype: str or :class:`~Fred2.Core.Protein.Protein`
        :raises ValueError: If start is greater than stop of index
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
            for pos, vs in sorted(self.vars.items()):
                if pos < start:
                    for v in vs:
                        if v.type in [VariationType.FSINS, VariationType.FSDEL]:
                            shift = (v.get_shift()+shift) % 3
                            if shift:
                                _fs.setdefault(pos-start, []).append(v)
                            else:
                                _fs.clear()
                    if pos in slice:
                        _vars[pos-start] = vs
            _vars.update(_fs)
            trans_id = self.transcript_id+":"+str(Protein.newid())
            seq = str(self)[index]
            return Protein(seq, gene_id=self.gene_id, transcript_id=trans_id, vars=_vars)

    def __repr__(self):
        # Header:
        lines = []
        lines += ["PROTEIN: %s (aa-seq)" % str(self)]
        lines += ["\t  %s (orig transcript)"%self.transcript_id]

        # Variants:
        lines += ["\t VARIANTS:"]
        for vpos, vset in self.vars.items():
            for v in vset:
                lines.append('\t pos %i: %s'%(vpos, v))

        return '\n\t'.join(lines) + '\n'

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

