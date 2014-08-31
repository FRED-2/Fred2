# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Peptide
   :synopsis: Contains the Peptide Class.
.. moduleauthor:: brachvogel

"""
__author__ = 'brachvogel'

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Fred2.Core.Base import MetadataLogger


class Peptide(MetadataLogger, Seq):
    """
    Peptides, belonging to one or several proteins.

    .. note::

        For accessing and manipulating the sequence see also :mod:`Bio.Seq.Seq`
        (from Biopython)

    :param dict(str,Protein) proteins: dict of transcript_IDs to protein
                                       instances that could generate that 
                                       peptide
    :param dict(int,Variant) vars: dict of positions to variant instances that
                                   affeced the peptide
    :param dict(str,Transcript) transcripts: dict of transcript_IDs to 
                                             transcript instances that could 
                                             have generated the peptide
    """

    def __init__(self, _seq):
        """
        :param str _seq: sequence of the peptide in one letter amino acid code
        """
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq, IUPAC.IUPACProtein)
        self.proteins = {}
        self.vars = dict()
        self.transcripts = {}


    def __getitem__(self, index):
        """
        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)
        
        :param int index: position in the peptide sequence
        :returns: (Peptide) -- A Peptide consisting of the single letter at
                  position :attr:`index`.
        """
        item = self[index]
        new_pept = Peptide(item)
        new_pept.proteins = self.proteins
        new_pept.vars = self.vars
        new_pept.transcripts = self.transcripts
        return new_pept


    def __repr__(self):
        lines = []
        for vpos, vset in self.vars.iteritems():
            lines.append('%s: '%vpos +', '.join([('%s %s' % \
            (v.coding.protPos, v.coding.aaMutationSyntax)) for v in vset]))
        lines.append("found in transcripts: \n" + ', '.join(\
            [t_ids for t_ids in self.transcripts]))
        return str(self) + '\n\t' + '\n\t'.join(lines)

