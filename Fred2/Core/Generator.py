# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.Generator
   :synopsis: Contains functions to transform variants to transcripts and 
              proteins to peptides. Transition of transcripts to proteins
              is done via :meth:`~Fred2.Core.Transcript.translate`
.. moduleauthor:: schubert

"""

import warnings
import collections

from Fred2.Core.Base import COMPLEMENT
from Fred2.Core.Protein import Protein
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Transcript import Transcript
from Fred2.Core.Variant import VariationType
from Fred2.IO.ADBAdapter import ADBAdapter, EAdapterFields

################################################################################
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

    :param list(char) seq: transcript sequence as a list
    :param Variant var: the snp variant to incorporate
    :param str transId: the transcript ID of seq
    :param int offset: the offset which has to be added onto the transcript
                       position of the variant
    :return: (list,int) the modified seq, the modified offset
    """
    if VariationType.SNP != var.type:
        raise TypeError("%s is not a SNP"%str(var))
    var.offsets[transId] = offset

    #print transId, len(seq), var.get_transcript_position(transId)-1
    if seq[var.get_transcript_position(transId)-1] != var.ref:
        warnings.warn("For %s bp dos not mmatch ref of assigned variant %s. Pos %i, var ref %s, seq ref %s " % (
        transId, str(var), var.get_transcript_position(transId) - 1, var.ref,
        seq[var.get_transcript_position(transId) - 1]))

    seq[var.get_transcript_position(transId)-1] = var.obs

    return seq, offset


def _incorp_insertion(seq, var, transId, offset):
    """
    incorporates an insertion into the given transcript sequence

    :param list(char) seq: (list) transcript sequence as a list
    :param Variant var: the snp variant to incorporate
    :param str transId: the transcript ID of seq
    :param int offset: the offset which has to be added onto the transcript
                       position of the variant
    :return: (list,int) modified sequence, the modified offset
    """
    if var.type not in [VariationType.INS, VariationType.FSINS]:
        raise TypeError("%s is not a insertion"%str(var))

    var.offsets[transId] = offset
    pos = var.get_transcript_position(transId)

    seq[pos:pos] = var.obs
    return seq, offset + len(var.obs)


def _incorp_deletion(seq, var, transId, offset):
    """
    incorporates a deletion into the given transcript sequence

    :param list(char) seq: transcript sequence as a list
    :param Variant var: the snp variant to incorporate
    :param str transId: the transcript ID of seq
    :param int offset: the offset which has to be added onto the transcript 
                       position of the variant
    :return: (list,int) -- modified sequence, the modified offset
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


def _check_for_problematic_variants(vars):
    """
    Filters problematic variants, e.g. variants that coincide.

    :param list(Variant) vars: initial list of variants
    :return: boole -- ture if now intersecting variants were found
    :invariant: list(Variant) vars: List is sorted based on genome position
    """
    v_tmp = vars[:]
    v = v_tmp.pop()
    current_range = (v.genomePos, v.genomePos
                                      +len(v.ref)-1 if v.type in [VariationType.FSDEL, VariationType.DEL] else
                                      v.genomePos)
    for v in reversed(v_tmp):
        if v.genomePos <= current_range[1]:
            #print "crash",current_range, v
            return False
        else:
            current_range = (v.genomePos, v.genomePos
                                      +len(v.ref)-1 if v.type in [VariationType.FSDEL, VariationType.DEL] else
                                      v.genomePos)
            #print "new block",v, current_range
    return True


#################################################################
# Public transcript generator functions

################################################################################
#        V A R I A N T S     = = >    P E P T I D E S
################################################################################
def generate_peptides_from_variants(vars, length, dbadapter):
    """
    generates polymorphic peptides based on given variants

    :param list(Variant) vars: list of Variants
    :param int length: length of peptides
    :param dbadapter: DBAdapter to fetch the transcript sequences
    :return: List(Peptide) -- A list of polymorphic peptides
    """
    def _generate_combinations(tId, vs, seq, usedVs, offset):
        """
        recursive variant combination generator
        """
        transOff = generate_peptides_from_variants.transOff
        #print "TransOffset ", transOff, tId,usedVs
        if vs:
            v = vs.pop()
            if v.isHomozygous:
                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)
                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs[:]

                for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset):
                    yield s

                # update the transcript variant id
                old_trans = generate_peptides_from_variants.transOff
                generate_peptides_from_variants.transOff += 1
                transOff = generate_peptides_from_variants.transOff
                _update_var_offset(usedVs, tId+":FRED2_%i"%old_trans, tId+":FRED2_%i"%(transOff))

                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)

                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
        else:
            yield tId+":FRED2_%i"%transOff, seq, usedVs

    def _generate_heterozygous(tId, vs, seq, usedVs, offset=None):
        """
            incorporates heterozygous variants into polymorphic sequences
        """
        if vs:
            v = vs.pop()
            #find offset
            offset = 0
            if offset is None:
                for uv in usedVs:
                    if uv.genomePos < v.genomePos:
                        offset = uv.offsets[tId]
                    else:
                        break

            #generate combinatorial branches:
            vs_tmp = vs[:]
            tmp_seq = seq[:]
            tmp_usedVs = usedVs[:]

            for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset=offset):
                    yield s

            old_trans = generate_peptides_from_variants.transOff
            generate_peptides_from_variants.transOff += 1
            transOff = generate_peptides_from_variants.transOff
            _update_var_offset(usedVs, tId+":FRED2_%i"%old_trans, tId+":FRED2_%i"%(transOff))

            seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)

            usedVs.append(v)
            for s in _generate_combinations(tId, vs, seq, usedVs, offset=offset):
                yield s
        else:
            yield tId+":FRED2_%i"%generate_peptides_from_variants.transOff, seq, usedVs

    if not isinstance(dbadapter, ADBAdapter):
        raise ValueError("The given dbadapter is not of type ADBAdapter")

    transToVar = {}
    for v in vars:
        for trans_id in v.coding.iterkeys():
            transToVar.setdefault(trans_id, []).append(v)

    for tId, vs in transToVar.iteritems():
        print tId
        query = dbadapter.get_transcript_information(tId)
        if query is None:
            warnings.warn("Transcript with ID %s not found in DB"%tId)
            continue

        tSeq = query[EAdapterFields.SEQ]
        geneid = query[EAdapterFields.GENE]
        strand = query[EAdapterFields.STRAND]


        #if its a reverse transcript form the complement of the variants
        if strand == "-":
            for v in transToVar[tId]:
                v.ref = v.ref[::-1].translate(COMPLEMENT)
                v.obs = v.obs[::-1].translate(COMPLEMENT)

        vs.sort(key=lambda v: v.genomePos, reverse=True)
        if not _check_for_problematic_variants(vs):
            warnings.warn("Intersecting variants found for Transcript %s"%tId)
            continue
        generate_peptides_from_variants.transOff = 0
        prots = []
        vs_homo_and_fs = filter(lambda x: x.type in [VariationType.FSINS, VariationType.FSDEL] or x.isHomozygous, vs)
        vs_hetero = filter(lambda x: not x.isHomozygous, vs)

        prots = []
        for tId, varSeq, varComb in _generate_combinations(tId, vs_homo_and_fs, list(tSeq), [], 0):
            if vs_hetero:
                for i in xrange(len(varSeq)+1-3*length):
                    end = i+3*length
                    frac_seq = varSeq[i:end]
                    frac_var = filter(lambda x: i <= x.get_transcript_position < end, vs_hetero)
                    for ttId, vvarSeq, vvarComb in _generate_heterozygous(tId, frac_var, frac_seq, varComb):
                        prots.append(Transcript("".join(vvarSeq), geneid, ttId, _vars=vvarComb).translate())
            else:
                prots.append(Transcript("".join(varSeq), geneid, tId, _vars=varComb).translate())

        return generate_peptides_from_protein(prots, length)

################################################################################
#        V A R I A N T S     = = >    T R A N S C R I P T S
################################################################################


def generate_transcripts_from_variants(vars, dbadapter):
    """
    generates all possible transcript variations of the given variants

    :param list(Variant) vars: A list of variants for which transcripts should 
                               be build
    :param ADBAdapter dbadapter: a DBAdapter to fetch the transcript sequences
    :return: (Generator(Transcripts)) -- a generator of transcripts with all 
             possible variations determined by the given
             variant list
    """
    def _generate_combinations(tId, vs, seq, usedVs, offset):
        """
        recursive variant combination generator
        """
        transOff = generate_transcripts_from_variants.transOff
        #print "TransOffset ", transOff, tId,usedVs
        if vs:
            v = vs.pop()
            if v.isHomozygous:
                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)
                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs[:]

                for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset):
                    yield s

                # update the transcript variant id
                old_trans = generate_transcripts_from_variants.transOff
                generate_transcripts_from_variants.transOff += 1
                transOff = generate_transcripts_from_variants.transOff
                _update_var_offset(usedVs, tId+":FRED2_%i"%old_trans, tId+":FRED2_%i"%(transOff))

                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)

                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
        else:
            yield tId+":FRED2_%i"%transOff, seq, usedVs

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
        print tId
        query = dbadapter.get_transcript_information(tId)
        if query is None:
            warnings.warn("Transcript with ID %s not found in DB"%tId)
            continue

        tSeq = query[EAdapterFields.SEQ]
        geneid = query[EAdapterFields.GENE]
        strand = query[EAdapterFields.STRAND]


        #if its a reverse transcript form the complement of the variants
        if strand == "-":
            for v in transToVar[tId]:
                v.ref = v.ref[::-1].translate(COMPLEMENT)
                v.obs = v.obs[::-1].translate(COMPLEMENT)

        vs.sort(key=lambda v: v.genomePos, reverse=True)
        if not _check_for_problematic_variants(vs):
            warnings.warn("Intersecting variants found for Transcript %s"%tId)
            continue
        generate_transcripts_from_variants.transOff = 0
        for tId, varSeq, varComb in _generate_combinations(tId, vs, list(tSeq), [], 0):
            yield Transcript("".join(varSeq), geneid, tId, _vars=varComb)


def generate_transcripts_from_tumor_variants(normal, tumor, dbadapter):
    """
    generates all possible transcript variations of the given variants

    :param list(Variant) normal: A list of variants of the normal tissue
    :param list(Variant) tumor: A list of variant of the cancer tissue for which transcript should be generated
    :param ADBAdapter dbadapter: a DBAdapter to fetch the transcript sequences
    :return: (Generator(Transcripts)) -- a generator of transcripts with all
             possible variations determined by the given
             variant list
    """
    def _generate_combinations(tId, vs, seq, usedVs, offset):
        """
        recursive variant combination generator
        """
        transOff = generate_transcripts_from_tumor_variants.transOff
        #print "TransOffset ", transOff, tId,usedVs
        if vs:
            flag, v = vs.pop()

            if v.isHomozygous:
                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)
                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs[:]

                if flag:
                    for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset):
                        yield s

                # update the transcript variant id
                old_trans = generate_transcripts_from_tumor_variants.transOff
                generate_transcripts_from_tumor_variants.transOff += 1
                transOff = generate_transcripts_from_tumor_variants.transOff
                _update_var_offset(usedVs, tId+":FRED2_%i"%old_trans, tId+":FRED2_%i"%(transOff))

                seq, offset = _incorp.get(v.type, lambda a, b, c, d: (a, d))(seq, v, tId+":FRED2_%i"%transOff, offset)

                usedVs.append(v)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
        else:
            yield tId+":FRED2_%i"%transOff, seq, usedVs

    #1) get all transcripts and sort the variants to transcripts

    #For a transcript do:
        #A) get transcript sequences
        #B) generate all possible combinations of variants
        #C) apply variants to transcript and generate transcript object

    if not isinstance(dbadapter, ADBAdapter):
        raise ValueError("The given dbadapter is not of type ADBAdapter")

    transToVar = {}
    for v in tumor:
        for trans_id in v.coding.iterkeys():
            transToVar.setdefault(trans_id, []).append((False, v))

    for v in normal:
        for trans_id in v.coding.iterkeys():
            if trans_id in transToVar:
                transToVar.setdefault(trans_id, []).append((True, v))

    for tId, vs in transToVar.iteritems():
        query = dbadapter.get_transcript_information(tId)
        if query is None:
            warnings.warn("Transcript with ID %s not found in DB"%tId)
            continue

        tSeq = query[EAdapterFields.SEQ]
        geneid = query[EAdapterFields.GENE]
        strand = query[EAdapterFields.STRAND]


        #if its a reverse transcript form the complement of the variants
        if strand == "-":
            for flag, v in transToVar[tId]:
                v.ref = v.ref[::-1].translate(COMPLEMENT)
                v.obs = v.obs[::-1].translate(COMPLEMENT)

        vs.sort(key=lambda v: v[1].genomePos, reverse=True)
        if not _check_for_problematic_variants(map(lambda x: x[1],vs)):
            warnings.warn("Intersecting variants found for Transcript %s"%tId)
            continue

        generate_transcripts_from_tumor_variants.transOff = 0
        for tId, varSeq, varComb in _generate_combinations(tId, vs, list(tSeq), [], 0):
            yield Transcript("".join(varSeq), geneid, tId, _vars=varComb)


################################################################################
#        P R O T E I N    = = >    P E P T I D E
################################################################################

def generate_peptides_from_protein(proteins, window_size):
    """
    Creates all peptides for a given window size, from a given protein. The
    result is a generator.

    :param Protein protein: (list of) protein(s) from which a list of unique
                            peptides should be generated
    :param int window_size: size of peptide fragments
    """
    def frameshift_influences(tid, _vars, res, start):
        # find variants out side the peptide frame, still influencing it via a
        # frameshift
        accu = [] # accumulator for relevant variants

        _vars.sort(key=lambda v: v.genomePos) # necessary?
        shift = 0

        for var in _vars:

            pos = var.get_protein_position(tid)
            new_shift = var.get_shift()

            if pos < start:
                # does a variant yield a frame shift?
                if shift + new_shift:
                    shift += new_shift
                    accu.append(var)
                else:
                    accu = {}
            # here: var.get_protein_position >= start, we are done!
            else:
                res += accu
                break

    def gen_peptide_info(protein):
        # Generate peptide sequences and find the variants within each
        res = []

        seq = str(protein)
        for i in xrange(len(protein)+1-window_size):
            # generate peptide fragment
            end = i+window_size
            pep_seq = seq[i:end]

             # get the variants affecting the peptide:
            if protein.vars:
                # variants within the peptide:
                pep_var = [var for pos, var_list in protein.vars.iteritems() \
                           for var in var_list if i <= pos <= end]

                # outside variants that affect the peptide via frameshift:
                frameshift_influences(protein.transcript_id, 
                                      protein.orig_transcript.vars.values(),
                                      pep_var, i)
            else:
                pep_var = []

            res.append((pep_seq, pep_var))
        return res

    final_peptides = {} # sequence : peptide-instance

    if isinstance(proteins, Protein):
        proteins = [proteins]

    if any(not isinstance(p, Protein) for p in proteins):
        raise ValueError("Input does contain non protein objects.")

    for prot in proteins:
        # generate all peptide sequences per protein:
        for (seq, _vars) in gen_peptide_info(prot):

            t_id = prot.transcript_id
            if seq not in final_peptides:
                final_peptides[seq] = Peptide(seq)

            final_peptides[seq].proteins[t_id] = prot
            final_peptides[seq].vars[t_id] = _vars
            final_peptides[seq].transcripts[t_id] = prot.orig_transcript

    return final_peptides.values()
