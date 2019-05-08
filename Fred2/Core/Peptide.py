# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.Peptide
   :synopsis: Contains the Peptide class

   :Note: All internal indices start at 0!

.. moduleauthor:: schubert,walzer

"""
import collections

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Fred2.Core import MetadataLogger
from Fred2.Core.Protein import Protein
from Fred2.Core.Variant import VariationType


class Peptide(MetadataLogger, Seq):
    """
    This class encapsulates a :class:`~Fred2.Core.Peptide.Peptide`, belonging to one or several
    :class:`~Fred2.Core.Protein.Protein`.

    .. note:: For accessing and manipulating the sequence see also :mod:`Bio.Seq.Seq` (from Biopython)
    """

    def __init__(self, seq, protein_pos=None):
        """
        :param str seq: Sequence of the peptide in one letter amino acid code
        :param protein_pos: Dict of transcript_IDs to position of origin in protein
        :type protein_pos: dict(:class:`~Fred2.Core.Protein.Protein`,list(int))`

        """
        MetadataLogger.__init__(self)
        Seq.__init__(self, seq.upper(), IUPAC.IUPACProtein)

        # Enforce dict storage
        if protein_pos and \
                any(not isinstance(p, Protein) or any(not isinstance(i, int) for i in pos) for p, pos in
                    protein_pos.items()):
            raise TypeError("The proteins_pos given to a Peptide object should be dict(Protein,list(int))")
        self.proteins = dict() if protein_pos is None else {p.transcript_id:p for p in protein_pos.keys()}
        self.proteinPos = collections.defaultdict(list) if protein_pos is None else {p.transcript_id: pos for p, pos in
                                                                                     protein_pos.items()}

    def __getitem__(self, index):
        """
        Overrides :meth:`Bio.Seq.Seq.__getitem__` (from Biopython)

        Returns a single letter or a sliced :class:`~Fred2.Core.Peptide.Peptide`.
        Allows only simple slicing (i.e. start < stop)

        :param int/Slice index: position in the peptide sequence
        :return: A single letter at position :attr:`index` or a sliced :class:`~Fred2.Core.Peptide.Peptide`
        :rtype: :class:`~Fred2.Core.Peptide.Peptide`
        :raises ValueError: If stop is greater than start of index
        """
        if isinstance(index, int):
            #Return a single letter as a string
            return str(self)[index]
        else:
            seq = str(self)[index]
            start, stop, step = index.indices(len(self))
            if start > stop:
                raise ValueError("start has to be greater than stop")
            protPos = {self.proteins[tId]: [p+(start-p) for p in pos] for tId, pos in self.proteinPos.items()}
            return Peptide(seq, protein_pos=protPos)

    def __repr__(self):

        lines = ["PEPTIDE:\n %s" % str(self)]
        #http://stackoverflow.com/questions/1436703/difference-between-str-and-repr-in-python/2626364#2626364
        for t in self.get_all_transcripts():
            if t is not None:
                t_id = t.transcript_id if t is not None else ""
                lines.append("in TRANSCRIPT: %s" % t_id)
                lines.append("\tVARIANTS:")
                for var in self.get_variants_by_protein(t_id):
                    lines.append("\t%s" % var)
        for p in self.proteins.values():
            p_id = p.transcript_id
            lines.append("in PROTEIN: %s" % p_id)
        return '\n'.join(lines)

    def get_all_proteins(self):
        """
        Returns all :class:`~Fred2.Core.Protein.Protein` objects associated with the
        :class:`~Fred2.Core.Peptide.Peptide`

        :return: A list of :class:`~Fred2.Core.Protein.Protein`
        :rtype: list(:class:`~Fred2.Core.Protein.Protein`)
        """
        return list(self.proteins.values())

    def get_protein(self, transcript_id):
        """
        Returns a specific protein object identified by a unique transcript-ID

        :param str transcript_id: A :class:`~Fred2.Core.Transcript.Transcript` ID
        :return: A :class:`~Fred2.Core.Protein.Protein`
        :rtype: :class:`~Fred2.Core.Protein.Protein`
        """
        #Default return value None if peptide does not origin from protein?
        return self.proteins.get(transcript_id, None)

    def get_all_transcripts(self):
        """
        Returns a list of :class:`~Fred2.Core.Transcript.Transcript` objects that are associated with the
        :class:`~Fred2.Core.Peptide.Peptide`

        :return: A list of :class:`~Fred2.Core.Transcript.Transcript`
        :rtype: list(:class:`~Fred2.Core.Transcript.Transcript`)
        """
        return [p.orig_transcript for p in self.proteins.values()]

    def get_transcript(self, transcript_id):
        """
        Returns a specific :class:`~Fred2.Core.Transcript.Transcript` object identified by a unique transcript-ID

        :param str transcript_id: A :class:`~Fred2.Core.Transcript.Transcript` ID
        :return: A :class:`~Fred2.Core.Transcript.Transcript`
        :rtype: :class:`~Fred2.Core.Transcript.Transcript`
        """
        try:
            return self.proteins[transcript_id].orig_transcript
        except KeyError:
            return None

    def get_protein_positions(self, transcript_id):
        """
        Returns all positions of origin for a given :class:`~Fred2.Core.Protein.Protein` identified by its transcript-ID

        :param str transcript_id: The unique transcript ID of the :class:`~Fred2.Core.Protein.Protein` in question
        :return: A list of positions within the protein from which the :class:`~Fred2.Core.Peptide.Peptide`
                 originated (starts at 0)
        :rtype:  list(int)
        """
        return self.proteinPos.get(transcript_id, [])

    def get_variants_by_protein(self, transcript_id):
        """
        Returns all :class:`~Fred2.Core.Variant.Variant` of a :class:`~Fred2.Core.Protein.Protein` that have influenced
        the :class:`~Fred2.Core.Peptide.Peptide` sequence

        :param str transcript_id: :class:`~Fred2.Core.Transcript.Transcript` ID of the specific protein in question
        :return: A list variants that influenced the peptide sequence
        :rtype: list(:class:`~Fred2.Core.Variant.Variant`)
        :raises KeyError: If peptide does not originate from specified :class:`~Fred2.Core.Protein.Protein`
        """
        try:
            p = self.proteins[transcript_id]
            var = []
            fs = []
            shift = 0
            for start_pos in self.proteinPos[transcript_id]:
                for i in range(start_pos):
                    for v in p.vars.get(i, []):
                        if v.type in [VariationType.FSDEL, VariationType.FSINS]:
                            shift = (v.get_shift()+shift) % 3
                            if shift:
                                fs.append(v)
                            else:
                                fs = []
                for j in range(start_pos, start_pos+len(self)):
                    for v in p.vars.get(j, []):
                        var.append(v)
            fs.extend(var)
            return fs
        except KeyError:
            raise KeyError("Peptide does not origin from protein with \
            transcript ID {transcript}".format(transcript=transcript_id))

    def get_variants_by_protein_position(self, transcript_id, protein_pos):
        """
        Returns all :class:`~Fred2.Core.Variant.Variant` and their relative position to the peptide sequence of a given
        :class:`~Fred2.Core.Protein.Protein` and protein position

        :param str transcript_id: A :class:`~Fred2.Core.Transcript.Transcript` ID of the specific protein in question
        :param int protein_pos: The :class:`~Fred2.Core.Protein.Protein` position at which the peptides sequence starts
                                in the protein
        :return: Dictionary of relative position of variants in peptide (starts at 0) and associated variants that
                 influenced the peptide sequence
        :rtype: dict(int,list(:class:`~Fred2.Core.Variant.Variant`))
        :raises:
         :ValueError: If :class:`~Fred2.Core.Peptide.Peptide` does not start at specified position
         :KeyError: If :class:`~Fred2.Core.Peptide.Peptide` does not originate from specified
                    :class:`~Fred2.Core.Protein.Protein`
        """
        try:
            p = self.proteins[transcript_id]
            if protein_pos not in self.proteinPos[transcript_id]:
                    raise ValueError("Peptide does not start a \
                                     {pos} in protein with transcript ID {transcript}".format(pos=protein_pos,
                                                                                              transcript=protein_pos))
            var = dict()
            fs = dict()
            shift = 0
            for i in range(protein_pos):
                for v in p.vars.get(i, []):
                    if v.type in [VariationType.FSDEL, VariationType.FSINS]:
                        shift = (v.get_shift()+shift) % 3
                        if shift:
                            fs.setdefault(i - protein_pos, []).append(v)
                        else:
                            fs.clear()
            for j in range(protein_pos, protein_pos+len(self)):
                for v in p.vars.get(j, []):
                    var.setdefault(j - protein_pos, []).append(v)
            fs.update(var)
            return fs
        except KeyError:
            raise KeyError("Peptide does not origin from protein with \
                             transcript ID {transcript}".format(transcript=transcript_id))

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
