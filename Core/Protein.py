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

from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import MetadataLogger




class Protein(MetadataLogger, Seq):

    def __init__(self, _seq, _gene_id, _transcript_id, _orig_transcript=None, 
                 _vars=None):
        """
        :param _seq: String of an IUPACProtein alphabet, representing a protein
        :type _sep: str.
        :param _gene_id: ID of the genome
        :type _gene_id: str.
        :param _orig_transcript: Reference to the originating transcript
        :type _orig_transcript: Fred2.Core.Transcript.Transcript
        :param _vars: Nonsynonymous(!) variants that are assoziated with the
                      protein/transcript. Key: their protein position
        :type _vars: {int : Fred2.Core.Variant.Variant}
        """
        # Init parent type:
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq, IUPAC.IUPACProtein)
        # Init own member:
        if _vars is None:
            self.vars = dict()
        else:
            self.vars = _vars  # {prot-position: variantset}
        self.orig_transcript = _orig_transcript
        self.transcript_id = _transcript_id
        self.gene_id = _gene_id


    def __getitem__(self, index):
        """
        :param index: (int) position within the primary sequence
        :returns: Protein -- A protein consisting of the single letter at
            position :attr:`index`.
        """
        item = self[index]
        return Protein(item, self.gene_id, self.orig_transcript, self.vars)

    def __repr__(self):
        lines = ["sequ:"+str(self)] # the prot sequence
        for vpos, vset in self.vars.iteritems():
            lines.append('var at %s: '%vpos +', '.join([('%s %s' % \
            (v.coding.protPos, v.coding.aaMutationSyntax)) for v in vset]))
        return self.transcript_id + ', ' + ', '.join(lines)


def generate_peptides_from_protein(proteins, window_size):
    """
    Creates all peptides for a given window size, from a given protein. The
    result is a generator.
    :param protein:
    :type protein: Fred2.Core.Protein.Protein.
    :param window_size: size of peptide fragments
    :type window_size: int
    """

    def gen_peptide_info(protein):
        # Generate a peptide seqs and find the variants within each sequence

        res = []
        seq = str(protein)
        for i in range(len(protein)+1-window_size):
            # generate peptide fragment
            end = i+window_size
            pep_seq = seq[i:end]

             # get the variants within that peptide:
            if protein.vars is not None:
                pep_var = dict((pos, var) for pos, var in \
                            protein.vars.iteritems() if i <= pos <= end)
            else:
                pep_var = None

            res.append((pep_seq, pep_var))
        return res


    final_peptides = {} # sequence : peptide-instance

    for prot in proteins:
        # generate all peptide sequences per protein:
        for (pep, _vars) in gen_peptide_info(prot):

            t_id = prot.transcript_id
            if pep not in final_peptides:
                final_peptides[pep] = Peptide(pep)

            final_peptides[pep].proteins[t_id] = prot
            final_peptides[pep].vars[t_id] = _vars
            final_peptides[pep].transcripts[t_id] = prot.orig_transcript

    return final_peptides.values()
