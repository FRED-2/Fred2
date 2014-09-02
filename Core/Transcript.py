# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Transcript
   :synopsis: Contains the Transcript Class.
.. moduleauthor:: brachvogel, szolek, walzer

"""
__author__ = 'brachvogel', 'szolek', 'walzer'

import warnings

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

    def __init__(self, _gene_id, _transcript_id, _seq, _vars=None):
        """
        :param str _gene_id: input genome ID
        :param str _transcript_id: input transcript RefSeqID
        :param str _seq: Transcript RefSeq sequence
        :param dict(int,Variant) _vars: Dict of Variants for specific positions 
                                        in the transcript. key=position, 
                                        value=Variant
        """
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq, generic_rna)
        self.gene_id = _gene_id
        self.transcript_id = _transcript_id
        if _vars is not None:
            self.vars = {v.get_transcript_position(_transcript_id): \
            v for v in _vars}
        else:
            self.vars = dict()


    def __getitem__(self, index):
        """
        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)

        :param int index: position 
        :returns: (Transcript) -- A Transcript consisting of the single
        letter at position :attr:`index`.
        """
        item = self[index]
        return Transcript(self.transcript_id, item, self.vars)


    def __repr__(self):
        # get sequence:
        lines = [str(self)]
        # get all variants:
        t_id = self.transcript_id.split(":FRED2_")[0]
        for vpos, vset in self.vars.iteritems():
            lines.append('%s: '%vpos +', ' + ('%s %s' % \
            (vset.get_transcript_position(self.transcript_id), \
            vset.coding[t_id].aaMutationSyntax)))

        return self.transcript_id + '\n\t' + '\n\t'.join(lines)


    def translate(self, table='Standard', stop_symbol='*', to_stop=False, 
                  cds=False):
        """
        Overrides :meth:`Bio.Seq.Seq.translate` (from Biopython) and enables 
        the translation from a transcript to a protein instance

        :param returns: (Protein) -- the protein that corresponds to the 
                        transcript
        """
        # translate to a protein sequence
        if len(str(self)) % 3 != 0:
            raise ValueError('ERROR while translating: lenght of transcript %s \
is no multiple of 3, the transcript is:\n %s' % (self.transcript_id, self))

        prot_seq = str(Seq.translate(self))

        # only transfer the non-synonymous variants to the protein as an
        # ordered dict, also translate into protein positions
        new_vars = dict()
        for var in self.vars.values():
            if not var.isSynonymous:
                pos = var.get_protein_position(self.transcript_id)
                new_vars.setdefault(pos, []).append(var)

        gene_id = self.gene_id
        return Protein(prot_seq, gene_id, self.transcript_id, self, new_vars)

