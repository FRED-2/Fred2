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
    :param dict(int,list(Variant)) vars: dict of positions to a list of variants
                                         that affeced the peptide, (including
                                         frameshifts that started not directly
                                         within the peptide)
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
        lines = ["peptide: %s"%str(self)]

        # print self.vars
        # for vpos, vset in self.vars.iteritems():
        #     lines.append("pos %i: %s"% vpos, vset)

        for t_id in self.transcripts:
            lines.append("transcript: %s"%t_id)
            lines.append("\t variants:")

            for (pos, _vars) in self.vars[t_id].items():
                lines.append("\t pos %d: %s" % (pos, ', '.join([str(var) for var in _vars])))

        return '\n\t'.join(lines) + '\n'

