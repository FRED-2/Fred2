# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Protein
   :synopsis: Contains the Protein Class.
.. moduleauthor:: brachvogel

"""
__author__ = 'brachvogel'

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Fred2.Core.Base import MetadataLogger




class Protein(MetadataLogger, Seq):
    """
    Protein corresponding to exactly one transcript.

    .. note::

        For accessing and manipulating the sequence see also :mod:`Bio.Seq.Seq`
        (from Biopython)

    :param str gene_id: ID of the genome
    :param str transcript_id: ID of the corresponding transcript 
    :param Transcript orig_transcript: Reference to the originating transcript
    :param dict(int,list(Variant)) _vars: Nonsynonymous variants that are
                                          assoziated with the protein. 
                                          key=position within protein, 
                                          value=list of variants at that pos
    """

    def __init__(self, _seq, _gene_id, _transcript_id, _orig_transcript=None, 
                 _vars=None):
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
        Seq.__init__(self, _seq, IUPAC.IUPACProtein)
        # Init own member:
        if _vars is None:
            self.vars = dict()
        else:
            self.vars = _vars  # {prot-position: list(variant)}
        self.orig_transcript = _orig_transcript
        self.transcript_id = _transcript_id
        self.gene_id = _gene_id


    def __getitem__(self, index):
        """
        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)

        :param int index: position within the primary sequence
        :returns: Protein -- A protein consisting of the single letter at
                  position :attr:`index`.
        """
        item = str(self)[index]
        return Protein(item, self.gene_id, self.orig_transcript, self.vars)

    def __repr__(self):
        lines = ["sequ:"+str(self)] # the prot sequence
        for vpos, vset in self.vars.iteritems():
            lines.append('var at %s: '%vpos +', '.join([('%s %s' % \
            (v.coding.protPos, v.coding.aaMutationSyntax)) for v in vset]))
        return self.transcript_id + ', ' + ', '.join(lines)


