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
    :param dict(str,list(Variant)) vars: dict of transcript_IDs to a list of
                                         variants that affected the peptide,
                                         (including frame shifts that started not
                                         directly within the peptide)
    :param dict(str,Transcript) transcripts: dict of transcript_IDs to 
                                             transcript instances that could 
                                             have generated the peptide
    """

    def __init__(self, _seq, proteins=None, vars=None, vars_position=None, transcripts=None):
        """
        :param str _seq: sequence of the peptide in one letter amino acid code
        """
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq, IUPAC.IUPACProtein)
        self.proteins = {} if proteins is None else proteins
        self.vars = {} if vars is None else vars
        self.transcripts = {} if transcripts is None else transcripts


    def __getitem__(self, index):
        """
        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)
        
        :param int index: position in the peptide sequence
        :returns: (Peptide) -- A Peptide consisting of the single letter at
                  position :attr:`index`.
        """
        item = str(self)[index]
        new_pept = Peptide(item)
        new_pept.proteins = self.proteins
        new_pept.vars = self.vars
        new_pept.transcripts = self.transcripts
        return new_pept


    def __repr__(self):
        lines = ["PEPTIDE: %s"%str(self)]

        for t_id in self.transcripts:
            lines.append("TRANSCRIPT: %s"%t_id)
            lines.append("\t VARIANTS:")

            for var in self.vars[t_id]:
                lines.append("\t %s" %var)

        return '\n\t'.join(lines) + '\n'

    def get_all_variants(self):
        return [var for var_list in self.vars.values() for var in var_list]

    def get_all_proteins(self):
        return self.proteins.values()

    def get_all_transcripts(self):
        return self.transcripts.values()

    def __eq__(self, other):
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))