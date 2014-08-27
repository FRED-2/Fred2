# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'schubert'

from Core.Variant import VariationType
from IO.ADBAdapter import ADBAdapter
#Private module functions. It should not be possible to import these!

def _incorp_snp(seq, var, transId, offset):
    """
    incorporates a snp into the given transcript sequence
    :param seq: (list) transcript sequence as a list
    :param var: (Variant) the snp variant to incorporate
    :param transId: (str) the transcript ID of seq
    :param offset: (int) the offset which has to be added onto the transcript position of the variant
    :return: (list) the modified seq, (int) the modified offset
    """
    if VariationType.SNP != var.type:
        raise TypeError("%s is not a SNP"%str(var))

    seq[var.get_transcript_position(transId)+offset] = var.ref
    return seq, offset


def _incorp_insertion(seq, var, transId, offset):
    """

    !!!Danger site-effects for sequence!

    incorporates a snp into the given transcript sequence
    :param seq: (list) transcript sequence as a list
    :param var: (Variant) the snp variant to incorporate
    :param transId: (str) the transcript ID of seq
    :param offset: (int) the offset which has to be added onto the transcript position of the variant
    :return: (list) modified sequence, (int) the modified offset
    """
    if var.type not in [VariationType.INS, VariationType.FSINS]:
        raise TypeError("%s is not a insertion"%str(var))

    pos = var.get_transcript_position(transId)
    seq[pos+1+offset:pos+1+offset] = var.obs
    return seq, offset + len(var.observed)


def _incorp_deletion(seq, var, transId, offset):
    """

    !!!Danger site-effects for sequence!

    incorporates a snp into the given transcript sequence
    :param seq: (list) transcript sequence as a list
    :param var: (Variant) the snp variant to incorporate
    :param transId: (str) the transcript ID of seq
    :param offset: (int) the offset which has to be added onto the transcript position of the variant
    :return: (list) modified sequence, (int) the modified offset
    """
    if var.type not in [VariationType.DEL, VariationType.FSDEL]:
        raise TypeError("%s is not a deletion"%str(var))

    pos = var.get_transcript_position(transId)
    s = slice(pos + offset, pos+len(var.ref) + offset)
    del seq[s]
    return seq, offset - len(var.ref)

_incorp={VariationType.DEL: _incorp_deletion,
         VariationType.FSDEL: _incorp_deletion,
         VariationType.INS: _incorp_insertion,
         VariationType.FSINS: _incorp_insertion,
         VariationType.SNP: _incorp_snp
         }

#################################################################
# Public transcript generator functions





def generate_transcripts_from_variants(vars, dbadapter):
    """
    generates all possible transcript variations of the given variants

    :param vars: (list(Variation)) A list of variants for which transcripts should be build
    :param dbadapter: (ADBAdapter) a DBAdapter to fetch the transcript sequences
    :return: (list(Transcripts)) a list of transcripts with all possible variations determined by the given variant list
    """
    def _generate_combinations(tId, vs, seq, offset):
        """
         recursive variant combination generator
        :param tId:
        :param vs:
        :param seq:
        :param offset:
        :return:
        """
        if not vs:
            yield seq

        v = vs.pop()
        tmp_seq, tmp_offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId, offset)
        if v.isHomozygous:
            for s in _generate_combinations(vs, tmp_seq, offset):
                yield s
        else:
            for s in _generate_combinations(vs, tmp_seq, tmp_offset):
                yield s
            for s in _generate_combinations(vs, seq, offset):
                yield s

    #1) get all transcripts and sort the variants to transcripts

    #For a transcript do:
        #A) get transcript sequences
        #B) generate all possible combinations of variants
        #C) apply variants to transcript and generate transcript object

    if not isinstance(ADBAdapter, dbadapter):
        raise ValueError("The given dbadapter is not of type ADBAdapter")

    transReturn = []

    transToVar = {}
    for v in vars:
        for trans_id in v.coding.iterkeys():
            transToVar.setdefault(trans_id, []).append(v)

    for tId, vs in transToVar.iteritems():
        geneName, tSeq = dbadapter.get_transcript_sequence()[0].items()

        if tSeq is None:
            raise KeyError("Transcript with ID %s not found in DB"%tId)

        vs = sorted(vs, key=lambda v: v.genome_pos)
        valid_trans = list(_generate_combinations(tId, vs, list(tSeq), 0))


########