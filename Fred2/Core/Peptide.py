# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

__author__ = 'brachvogel,walzer'

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import Fred2
from Fred2.Core.Base import MetadataLogger


class Peptide(MetadataLogger, Seq):
    """
    This class encapsulates a Peptide, belonging to one or several Proteins.

    .. note:: For accessing and manipulating the sequence see also :mod:`Bio.Seq.Seq` (from Biopython)
    """

    def __init__(self, _seq, proteins=None, variants=None, transcripts=None):
        """
        :param str _seq: sequence of the peptide in one letter amino acid code
        :param dict(str,Protein) proteins: dict of transcript_IDs to protein
                                       instances that could generate that
                                       peptide
        :param dict(str,list(Variant)) variants: dict of transcript_IDs to a list of
                                         variants that affected the peptide,
                                         (including frame shifts that started not
                                         directly within the peptide)
        :param dict(str,Transcript) transcripts: dict of transcript_IDs to
                                             transcript instances that could
                                             have generated the peptide
        """
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq, IUPAC.IUPACProtein)

        self.proteins = dict() if proteins is None else proteins
        # Enforce dict storage
        if proteins and not isinstance(proteins, dict):
            raise TypeError("The proteins given to a Peptide object should be dict(str,Protein)")
        if proteins and not all(isinstance(v, Fred2.Core.Protein) and isinstance(k, str) for k,v in proteins.iteritems()):
            raise TypeError("The proteins given to a Peptide object should be dict(str,Protein)")

        self.transcripts = dict() if transcripts is None else transcripts
        # Enforce dict storage
        if transcripts and not isinstance(transcripts, dict):
            raise TypeError("The transcripts given to a Peptide object should be dict(str,Transcript)")
        if transcripts and not all(isinstance(v, Fred2.Core.Transcript) and isinstance(k, str) for k,v in transcripts.iteritems()):
            raise TypeError("The proteins given to a Peptide object should be dict(str,Transcript)")

        self.variants = dict() if variants is None else variants
        # Enforce dict list storage
        if variants and not isinstance(variants, dict):
            raise TypeError("The variants given to a Peptide object should be dict(str,list(Variant))")
        if variants and not all(isinstance(v, list) and isinstance(k, str) for k,v in variants.iteritems()):
            raise TypeError("The variants given to a Peptide object should be dict(str,list(Variant))")
        if variants and not all(isinstance(var, Fred2.Core.Variant) for var_list in variants.values() for var in var_list):
            raise TypeError("The variants given to a Peptide object should be dict(str,list(Variant))")
        #TODO necessary to sanity check if all variants are registered with transcripts?
        #TODO move transcripts in front of variants in __init__ signature? make clear: transcripts before variants
        #TODO do we actually need variants in __init__ if variants are supposed to be registered in the transcripts?
        #TODO in the end OOP is about encapsulating data in a intelligible way - variant is good, keep transcript and protein 'loose' and drop most of the dict stuff??

        #TODO register dummy transcripts of the Variants if any not registered yet - otherwise repr() will break!!!
        for t_id in self.variants:
            if t_id not in self.transcripts:
                self.transcripts[t_id] = Fred2.Core.Transcript(_seq="", _transcript_id=t_id)


    def __getitem__(self, index):
        """
        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)
        
        :param int index: position in the peptide sequence
        :returns: A Peptide consisting of the single letter at position :attr:`index`.
        :rtype: Peptide
        """
        item = str(self)[index]
        new_pept = Peptide(item)
        new_pept.proteins = self.proteins
        new_pept.variants = self.variants
        new_pept.transcripts = self.transcripts
        return new_pept

    def __repr__(self):
        lines = ["PEPTIDE:\n %s" % str(self)]
        #http://stackoverflow.com/questions/1436703/difference-between-str-and-repr-in-python/2626364#2626364
        for t_id in self.transcripts:
            lines.append("in TRANSCRIPT: %s" % t_id)
            lines.append("\tVARIANTS:")
            for var in self.variants[t_id]:
                lines.append("\t%s" % var)
        for p_id in self.proteins:
            lines.append("in PROTEIN: %s" % p_id)
        return '\n'.join(lines)

    def get_all_variants(self):
        """
        :return: a concatenated list of all contained variants
        """
        return [var for var_list in self.variants.values() for var in var_list]

    def get_all_proteins(self):
        return self.proteins.values()

    def get_all_transcripts(self):
        return self.transcripts.values()

    def __eq__(self, other):
        return str(self) == str(other)

    def __cmp__(self, other):
        return cmp(str(self), str(other))

    def __hash__(self):
        return hash(str(self))
