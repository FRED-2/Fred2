# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'schubert'

from Fred2.Core.Transcript import Transcript
from Fred2.Core.Variant import VariationType
from Fred2.IO.ADBAdapter import ADBAdapter, EAdapterFields, COMPLEMENT
#Private module functions. It should not be possible to import these!

def _update_var_offset(vars, transId_old, transId_new):
    """
    :param var:
    :param transId_old:
    :param transId_new:
    """
    for v in vars:
        offset = v.offsets[transId_old]
        v.offsets[transId_new] = offset

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

    seq[var.get_transcript_position(transId)] = var.ref
    return seq, offset


def _incorp_insertion(seq, var, transId, offset):
    """
    incorporates a snp into the given transcript sequence
    :param seq: (list) transcript sequence as a list
    :param var: (Variant) the snp variant to incorporate
    :param transId: (str) the transcript ID of seq
    :param offset: (int) the offset which has to be added onto the transcript position of the variant
    :return: (list) modified sequence, (int) the modified offset
    """
    if var.type not in [VariationType.INS, VariationType.FSINS]:
        raise TypeError("%s is not a insertion"%str(var))

    var.offsets[transId] = offset
    pos = var.get_transcript_position(transId)
    seq[pos+1:pos+1] = var.obs
    return seq, offset + len(var.observed)


def _incorp_deletion(seq, var, transId, offset):
    """
    incorporates a snp into the given transcript sequence
    :param seq: (list) transcript sequence as a list
    :param var: (Variant) the snp variant to incorporate
    :param transId: (str) the transcript ID of seq
    :param offset: (int) the offset which has to be added onto the transcript position of the variant
    :return: (list) modified sequence, (int) the modified offset
    """
    if var.type not in [VariationType.DEL, VariationType.FSDEL]:
        raise TypeError("%s is not a deletion"%str(var))

    var.offsets[transId] = offset
    pos = var.get_transcript_position(transId)
    s = slice(pos, pos+len(var.ref))
    del seq[s]
    return seq, offset - len(var.ref)


_incorp = {
            VariationType.DEL: _incorp_deletion,
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
    :return: (Generator(Transcripts)) a generator of transcripts with all possible variations determined by the given
             variant list
    """
    def _generate_combinations(tId, vs, seq, usedVs, offset, transOff):
        """
         recursive variant combination generator
        :param tId:
        :param vs:
        :param seq:
        :param offset:
        :return:
        """
        if vs:
            v = vs.pop()

            if v.isHomozygous:
                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)
                usedVs.append(v)
                for s in _generate_combinations(vs, seq, usedVs, offset, transOff):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs[:]
                tmp_offset = offset

                for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, tmp_offset, transOff):
                    yield s

                _update_var_offset(usedVs, tId+":FRED2_%i"%transOff, tId+":FRED2_%i"%(transOff+1))
                transOff += 1
                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)
                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset, transOff):
                    yield s
        else:
            yield seq, usedVs

    #1) get all transcripts and sort the variants to transcripts

    #For a transcript do:
        #A) get transcript sequences
        #B) generate all possible combinations of variants
        #C) apply variants to transcript and generate transcript object

    if not isinstance(dbadapter, ADBAdapter):
        raise ValueError("The given dbadapter is not of type ADBAdapter")

    transToVar = {}
    for v in vars:
        for trans_id in v.coding.iterkeys():
            transToVar.setdefault(trans_id, []).append(v)

    for tId, vs in transToVar.iteritems():
        query = dbadapter.get_transcript_information(tId)
        tSeq = query[EAdapterFields.SEQ]
        geneid = query[EAdapterFields.GENE]
        strand = query[EAdapterFields.STRAND]

        #if its a reverse transcript form the complement of the variants
        if strand == "-":
            for v in transToVar[tId]:
                v.ref = v.ref[::-1].translate(COMPLEMENT)
                v.obs = v.obs[::-1].translate(COMPLEMENT)

        if tSeq is None:
            raise KeyError("Transcript with ID %s not found in DB"%tId)

        vs.sort(key=lambda v: v.genomePos, reverse=True)
        print vs
        for varSeq, varComb in _generate_combinations(tId, vs, list(tSeq), [], 0, 0):
            yield Transcript(geneid, tId, "".join(varSeq), _vars=varComb)
