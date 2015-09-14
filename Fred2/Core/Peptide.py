# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

__author__ = 'schubert,walzer'
import collections

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Fred2.Core import MetadataLogger
from Fred2.Core import Protein
from Fred2.Core.Variant import VariationType


class Peptide(MetadataLogger, Seq):
    """
    This class encapsulates a Peptide, belonging to one or several Proteins.

    .. note:: For accessing and manipulating the sequence see also :mod:`Bio.Seq.Seq` (from Biopython)
    """

    def __init__(self, _seq, protein_pos=None):
        """
        :param str _seq: sequence of the peptide in one letter amino acid code
        :param dict(Protein,list(int)) protein_pos: dict of transcript_IDs to position of origin in protein

        """
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq.upper(), IUPAC.IUPACProtein)

        # Enforce dict storage
        if protein_pos and any(not isinstance(p, Protein) or
                                any(not isinstance(i, int) for i in pos)
                               for p, pos in protein_pos.iteritems()):
            raise TypeError("The proteins_pos given to a Peptide object should be dict(Protein,list(int))")
        self.proteins = dict() if protein_pos is None else {p.transcript_id:p for p in protein_pos.iterkeys()}
        self.proteinPos = collections.defaultdict(list) if protein_pos is None else {p.transcript_id: pos for p, pos in
                                                                                     protein_pos.iteritems()}

    def __getitem__(self, index):
        #TODO: does not work that way! Reimplement!
        """
        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)

        Returns a single letter or a sliced Peptide with.
        Allows only simple slicing (i.e. start < stop)

        :param int/Slice index: position in the peptide sequence
        :returns: A single letter at position :attr:`index` or a sliced Peptide.
        :rtype: Peptide
        """
        if isinstance(index, int):
            #Return a single letter as a string
            return str(self)[index]
        else:
            seq = str(self)[index]
            start, stop, step = index.indices(len(self))
            if start > stop:
                raise ValueError("start has to be greater than stop")
            protPos = {self.proteins[tId]: [p+(start-p) for p in pos] for tId, pos in self.proteinPos.iteritems()}
            return Peptide(seq, protein_pos=protPos)

    def __repr__(self):

        lines = ["PEPTIDE:\n %s" % str(self)]
        #http://stackoverflow.com/questions/1436703/difference-between-str-and-repr-in-python/2626364#2626364
        for t in self.get_all_transcripts():
            t_id = t.transcript_id
            lines.append("in TRANSCRIPT: %s" % t_id)
            lines.append("\tVARIANTS:")
            for var in self.get_variants_by_protein(t_id):
                lines.append("\t%s" % var)
        for p in self.proteins:
            p_id = p.transcript_id
            lines.append("in PROTEIN: %s" % p_id)
        return '\n'.join(lines)

    def get_all_proteins(self):
        return self.proteins.values()

    def get_protein(self, _transcript_id):
        #Default return value None if peptide does not origin from protein?
        return self.proteins.get(_transcript_id, None)

    def get_all_transcripts(self):
        return [p.orig_transcript for p in self.proteins.itervalues()]

    def get_transcript(self, _transcript_id):
        try:
            return self.proteins[_transcript_id].orig_transcript
        except KeyError:
            return None

    def get_protein_positions(self, _transcript_id):
        """
        returns all positions of origin for a given protein
        identified by its transcript_id

        :param (str) _transcript_id: The unique transcript ID of the protein in question
        :return: list(int) - a list of positions within the protein from which the peptide originated (starts at 0)
        """
        return self.proteinPos.get(_transcript_id, [])

    def get_variants_by_protein(self, _transcript_id):
        """
        returns all variants of a protein that have influenced the peptide sequence

        :param _transcript_id: Transcript ID of the specific protein in question
        :return: list(Variant) - A list variants that influenced the peptide sequence
        """
        try:
            p = self.proteins[_transcript_id]
            var = []
            fs = []
            shift = 0
            for start_pos in self.proteinPos[_transcript_id]:
                for i in xrange(start_pos):
                    for v in p.vars.get(i, []):
                        if v.type in [VariationType.FSDEL, VariationType.FSINS]:
                            shift = (v.get_shift()+shift) % 3
                            if shift:
                                fs.append(v)
                            else:
                                fs = []
                for j in xrange(start_pos, start_pos+len(self)):
                    for v in p.vars.get(j, []):
                        var.append(v)
            return fs.extend(var)
        except KeyError:
            raise ValueError("Peptide does not origin from protein with "
                             "transcript ID {transcript}".format(transcript=_transcript_id))

    def get_variants_by_protein_position(self, _transcript_id, _protein_pos):
        """
        returns all variants and their relative position to the peptide sequence of a given protein
        and protein protein position

        :param _transcript_id:
        :param _protein_pos:
        :return: dict(int,list(Vars)) - dictionary of relative position of variants in peptide (starts at 0)
                                        and associated variants that influenced the peptide sequence
        """
        try:
            p = self.proteins[_transcript_id]
            if _protein_pos not in self.proteinPos[_transcript_id]:
                    raise ValueError("Peptide does not start a "
                                     "{pos} in protein with transcript ID {transcript}".format(pos=_protein_pos,
                                                                                               transcript=_protein_pos))
            var = dict()
            fs = dict()
            shift = 0
            for i in xrange(_protein_pos):
                for v in p.vars.get(i, []):
                    if v.type in [VariationType.FSDEL, VariationType.FSINS]:
                        shift = (v.get_shift()+shift) % 3
                        if shift:
                            fs.setdefault(i-_protein_pos, []).append(v)
                        else:
                            fs = dict()
            for j in xrange(_protein_pos, _protein_pos+len(self)):
                for v in p.vars.get(j, []):
                    var.setdefault(j, []).append(v)
            return fs.update(var)
        except KeyError:
            raise ValueError("Peptide does not origin from protein with "
                             "transcript ID {transcript}".format(transcript=_transcript_id))

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
