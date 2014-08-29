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

    def __init__(self, _seq):
        """
        :param _seq: String sequence in one letter amino acid code
        :type _seq: str.
        :param _proteins: dict of transcript_Ids to protein instances.
        :type _proteins: {str: Fred2.Core.Protein.Protein}.
        :param _vars: dictionary for variants according to position.
        :type _vars: {int:Fred2.Core.Variant.Variant}.
        :param _transcripts: dictionary for transcripts according to their id.
        :type _transcripts: {str: Fred2.Core.Transcript.Transcript}.
        """

        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq, IUPAC.IUPACProtein)
        self.proteins = {}
        self.vars = dict()
        self.transcripts = {}


    def __getitem__(self, index):
        """

        :param index: (int) position in the peptide sequence
        :returns: Peptide -- A Peptide consisting of the single letter at
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

