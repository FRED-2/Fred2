# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

__author__ = 'schubert,brachvogel,walzer,szolek'

import itertools

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Fred2.Core.Base import MetadataLogger
from Fred2.Core.Variant import VariationType


class Protein(MetadataLogger, Seq):
    """
    Protein corresponding to exactly one transcript.

    .. note::

        For accessing and manipulating the sequence see also :mod:`Bio.Seq.Seq`
        (from Biopython)

    :param str _gene_id: ID of the genome
    :param str _transcript_id: ID of the corresponding transcript
    :param Transcript _orig_transcript: Reference to the originating transcript
    :param dict(int,list(Variant)) _vars: Nonsynonymous variants that are
                                          assoziated with the protein. 
                                          key=position within protein, 
                                          value=list of variants at that pos
    """
    newid = itertools.count().next #this is evil and has no other purpose? it does not help that there may be more than one protein from one transcript - due to variants

    def __init__(self, _seq, _gene_id="unknown", _transcript_id=None, _orig_transcript=None, _vars=None):
        """
        :param str _seq: String of an IUPACProtein alphabet, representing the
                         protein
        :param str _gene_id: ID of the genome the protein originated from
        :param str _transcript_id: ID of the transcript the protein originated 
                                   from
        :param Transcript _orig_transcript: Reference to the originating 
                                            transcript
        :param dict(int,list(Variant)) _vars: Nonsynonymous variants that are
                                              assoziated with the protein. 
                                              key=position within protein, 
                                              value=list of variants at that pos
        """
        # Init parent type:
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq.upper(), IUPAC.IUPACProtein)
        # Init own member:
        if _vars is None:
            self.vars = dict()
        else:
            self.vars = _vars  # {prot-position: list(variant)}
        self.orig_transcript = _orig_transcript
        self.transcript_id = "Protein_%i"%Protein.newid() if _transcript_id is None else _transcript_id
        self.gene_id = _gene_id

    def __getitem__(self, index):
        """
        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)

        returns a single letter or a sliced protein (when given a slice).
        The sliced protein does not reference to a Transcript object!

        Allows only simple slicing (i.e. start < stop)

        :param int/Slice index: position within the primary sequence or a slice
        :returns: str/Protein - A single letter at position :attr:`index`
                  or a sliced Protein with adjusted variant positions.
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
            for pos, vs in sorted(self.vars.iteritems()):
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
            return Protein(seq, _gene_id=self.gene_id, _transcript_id=trans_id, _vars=_vars)

    def __repr__(self):
        # Header:
        lines = []
        lines += ["PROTEIN: %s (aa-seq)" % str(self)]
        lines += ["\t  %s (orig transcript)"%self.transcript_id]

        # Variants:
        lines += ["\t VARIANTS:"]
        for vpos, vset in self.vars.iteritems():
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

