# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Transcript
   :synopsis: Contains the Transcript Class.
.. moduleauthor:: brachvogel, szolek, walzer

"""
__author__ = 'brachvogel', 'szolek', 'walzer'

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Fred2.Core.Protein import Protein
from Fred2.Core.Base import MetadataLogger


class Transcript(MetadataLogger, Seq):
    """Transcript Class
    """

    def __init__(self, _gene_id, _transcript_id, _seq, _vars=None):
        """

        :param _transcript_id: Transcript RefSeqID
        :type _transcript_id: str.
        :param _seq: Transcript RefSeq sequence
        :type _seq: str.
        :param _start: true start position of transcript, positions starting at
                       1
        :type _start: int.
        :param _stop: true stop position of transcript
        :type _stop: int.
        :param _vars: List of variants for the transcript.
        :type _vars: [Fred2.Core.Variant].
        """
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq, IUPAC.IUPACUnambiguousRNA)
        self.gene_id = _gene_id
        self.transcript_id = _transcript_id
        if _vars is not None:
            self.vars = dict((v.get_transcript_position(_transcript_id), v) \
                for v in _vars)
        else:
            self.vars = dict()


    def __getitem__(self, index):
        """
        Overrides Bio.Seq.__getitem__

        :param index: (int) position 
        :returns: Transcript -- A Transcript consisting of the single letter at
            position :attr:`index`.
        """
        item = self[index]
        return Transcript(self.transcript_id, item, self.vars)


    def __repr__(self):
        lines = [str(self)]
        for vpos, vset in self.vars.iteritems():
            lines.append('%s: '%vpos +', '.join([('%s %s' % \
            (v.get_transcript_position(self.transcript_id), \
                v.coding[self.transcript_id].aaMutationSyntax)) for v in vset]))
        return self.transcript_id + '\n\t' + '\n\t'.join(lines)


    def translate(self, table='Standard', stop_symbol='*', to_stop=False, 
                  cds=False):
        """
        Override of Bio.Seq.translate()
        """

        # translate to a protein sequence
        prot_seq = str(self.translate())

        # only transfer the non-synonymous variants to the protein as an
        # ordered dict, also translate into protein positions
        new_vars = dict((y.coding.protPos, y) for y in \
                               self.vars.itervalues()  if y.isSynonymous)

        gene_id = self.gene_id
        return Protein(prot_seq, gene_id, self.transcript_id, self, new_vars)

