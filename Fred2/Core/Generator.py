# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.Generator
   :synopsis: Contains functions to transform variants to transcripts and 
              proteins to peptides.

   :Note: All internal indices start at 0!

.. moduleauthor:: schubert,walzer

"""

import warnings
from itertools import chain

from Fred2.Core.Base import COMPLEMENT
from Fred2.Core.Protein import Protein
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Transcript import Transcript
from Fred2.Core.Variant import VariationType
from Fred2.IO.ADBAdapter import ADBAdapter, EAdapterFields

################################################################################
#Private module functions. It should not be possible to import these!

#symbol for reverse complement
REVERS = "-"


def _incorp_snp(seq, var, transId, pos, offset):
    """
    incorporates a snp into the given transcript sequence

    :param list(char) seq: transcript sequence as a list
    :param Variant var: the snp variant to incorporate
    :param str transId: the transcript ID of seq
    :param int pos: the position of the variant
    :param int offset: the offset which has to be added onto the transcript
                       position of the variant
    :return: int - the the modified offset
    """
    if VariationType.SNP != var.type:
        raise TypeError("%s is not a SNP"%str(var))

    #print transId, len(seq), var.get_transcript_position(transId)-1
    if seq[pos] != var.ref:
        warnings.warn("For %s bp does not match ref of assigned variant %s. Pos %i, var ref %s, seq ref %s " % (
        transId, str(var), pos, var.ref,
        seq[pos]))

    seq[pos] = var.obs

    return offset


def _incorp_insertion(seq, var, transId, pos, offset):
    """
    incorporates an insertion into the given transcript sequence

    :param list(char) seq: (list) transcript sequence as a list
    :param Variant var: the snp variant to incorporate
    :param str transId: the transcript ID of seq
    :param int offset: the offset which has to be added onto the transcript
                       position of the variant
    :return: int - the modified offset
    """
    if var.type not in [VariationType.INS, VariationType.FSINS]:
        raise TypeError("%s is not a insertion"%str(var))

    seq[pos:pos] = var.obs
    return offset + len(var.obs)


def _incorp_deletion(seq, var, transId, pos, offset):
    """
    incorporates a deletion into the given transcript sequence

    :param list(char) seq: transcript sequence as a list
    :param Variant var: the snp variant to incorporate
    :param str transId: the transcript ID of seq
    :param int pos: the starting position of the deletion
    :param int offset: the offset which has to be added onto the transcript 
                       position of the variant
    :return: (list,int) -- modified sequence, the modified offset
    """
    if var.type not in [VariationType.DEL, VariationType.FSDEL]:
        raise TypeError("%s is not a deletion"%str(var))

    s = slice(pos, pos+len(var.ref))
    del seq[s]
    return offset - len(var.ref)


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
    :return: bool -- true if no intersecting variants were found
    :invariant: list(Variant) vars: List is sorted based on genome position in descending order
    """
    def get_range(var):
        current_range = [0,0]
        if var.type in [VariationType.FSDEL, VariationType.DEL]:
            current_range[0] = var.genomePos
            current_range[1] = var.genomePos+len(v.ref)-1
        elif var.type in [VariationType.FSINS, VariationType.INS]:
            current_range[0] = var.genomePos-1
            current_range[1] = var.genomePos-1
        else:
            current_range[0] = var.genomePos
            current_range[1] = var.genomePos
        return current_range
    v_tmp = vars[:]
    v = v_tmp.pop()

    current_range = get_range(v)
    for v in reversed(v_tmp):
        genome_pos = v.genomePos-1 if v.type in [VariationType.FSINS, VariationType.INS] else v.genomePos
        if genome_pos <= current_range[1]:
            return False
        else:
            current_range = get_range(v)
    return True


#################################################################
# Public transcript generator functions

################################################################################
#        V A R I A N T S     = = >    P E P T I D E S
################################################################################
def generate_peptides_from_variants(vars, length, dbadapter, peptides=None):
    """
    generates polymorphic peptides based on given variants

    :param list(Variant) vars: list of Variants
    :param int length: length of peptides
    :param dbadapter: DBAdapter to fetch the transcript sequences
    :param list(Peptide) peptides: a list of existing peptides. These will be updated with new variant information
                                   if a sequence contained in the list is generated by the new variants
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
                pos = v.coding[tId].tranPos + offset
                usedVs[pos] = v
                offset = _incorp.get(v.type, lambda a, b, c, d, e: e)(seq, v, tId, pos, offset)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs.copy()

                #generate transcript without the current variant
                for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset):
                    yield s

                #and one transcript with current variant as we can't resolve haplotypes
                # update the transcript variant id
                generate_peptides_from_variants.transOff += 1
                transOff = generate_peptides_from_variants.transOff
                pos = vs.coding[tId].tranPos + offset
                usedVs[pos] = v
                offset = _incorp.get(v.type, lambda a, b, c, d, e: e)(seq, v, tId, pos, offset)

                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
        else:
            yield tId + ":FRED2_%i"%transOff, seq, usedVs

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
            tmp_usedVs = usedVs.copy()

            for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset=offset):
                    yield s

            generate_peptides_from_variants.transOff += 1
            pos = v.coding[tId].tranPos + offset
            usedVs[pos] = v
            offset = _incorp.get(v.type, lambda a, b, c, d, e: e)(seq, v, tId, pos, offset)

            for s in _generate_combinations(tId, vs, seq, usedVs, offset=offset):
                yield s
        else:
            yield tId + ":FRED2_%i"%generate_peptides_from_variants.transOff, seq, usedVs

    if not isinstance(dbadapter, ADBAdapter):
        raise ValueError("The given dbadapter is not of type ADBAdapter")

    transToVar = {}
    for v in vars:
        for trans_id in v.coding.iterkeys():
            transToVar.setdefault(trans_id, []).append(v)

    for tId, vs in transToVar.iteritems():
        #print tId
        query = dbadapter.get_transcript_information(tId)
        if query is None:
            warnings.warn("Transcript with ID %s not found in DB"%tId)
            continue

        tSeq = query[EAdapterFields.SEQ]
        geneid = query[EAdapterFields.GENE]
        strand = query[EAdapterFields.STRAND]

        #if its a reverse transcript form the complement of the variants
        if strand == REVERS:
            for v in transToVar[tId]:
                v.ref = v.ref[::-1].translate(COMPLEMENT)
                v.obs = v.obs[::-1].translate(COMPLEMENT)

        vs.sort(key=lambda v: v.genomePos - 1
                if v.type in [VariationType.FSINS, VariationType.INS]
                else v.genomePos, reverse=True)
        if not _check_for_problematic_variants(vs):
            warnings.warn("Intersecting variants found for Transcript %s"%tId)
            continue
        generate_peptides_from_variants.transOff = 0
        vs_homo_and_fs = filter(lambda x: x.type in [VariationType.FSINS, VariationType.FSDEL] or x.isHomozygous, vs)
        vs_hetero = filter(lambda x: not x.isHomozygous
                           and x.type not in [VariationType.FSINS, VariationType.FSDEL], vs)

        prots = []
        for tId, varSeq, varComb in _generate_combinations(tId, vs_homo_and_fs, list(tSeq), {}, 0):
            if vs_hetero:
                for i in xrange(len(varSeq) + 1 - 3 * length):
                    end = i + 3 * length
                    frac_seq = varSeq[i:end]
                    trans_id = tId.split(":FRED2_")[0]
                    offset = sum(v.get_transcript_offset() for pos, v in varComb.iteritems() if i <= pos <= end)
                    frac_var = filter(lambda x: i <= x.coding[trans_id].transPos+offset < end, vs_hetero)
                    for ttId, vvarSeq, vvarComb in _generate_heterozygous(tId, frac_var, frac_seq, varComb):
                        chain(prots, generate_proteins_from_transcripts(Transcript("".join(vvarSeq), geneid, ttId,
                                                                                       _vars=vvarComb)))
            else:
                if not prots:
                    prots = generate_proteins_from_transcripts(
                        Transcript("".join(varSeq), geneid, tId, _vars=varComb))
                else:
                    chain(prots, generate_proteins_from_transcripts(
                        Transcript("".join(varSeq), geneid, tId, _vars=varComb)))

        return generate_peptides_from_protein(prots, length, peptides=peptides)

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
                pos = v.coding[tId].tranPos + offset
                usedVs[pos] = v
                offset = _incorp.get(v.type, lambda a, b, c, d, e: e)(seq, v, tId, pos, offset)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs.copy()

                for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset):
                    yield s

                # update the transcript variant id
                generate_transcripts_from_variants.transOff += 1
                pos = v.coding[tId].tranPos + offset
                usedVs[pos] = v

                offset = _incorp.get(v.type, lambda a, b, c, d, e: e)(seq, v, tId, pos, offset)

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
        #print tId
        query = dbadapter.get_transcript_information(tId)
        if query is None:
            warnings.warn("Transcript with ID %s not found in DB"%tId)
            continue

        tSeq = query[EAdapterFields.SEQ]
        geneid = query[EAdapterFields.GENE]
        strand = query[EAdapterFields.STRAND]

        #if its a reverse transcript form the complement of the variants
        if strand == REVERS:
            for v in transToVar[tId]:
                v.ref = v.ref[::-1].translate(COMPLEMENT)
                v.obs = v.obs[::-1].translate(COMPLEMENT)

        vs.sort(key=lambda v: v.genomePos-1
                if v.type in [VariationType.FSINS, VariationType.INS]
                else v.genomePos, reverse=True)
        if not _check_for_problematic_variants(vs):
            warnings.warn("Intersecting variants found for Transcript %s"%tId)
            continue
        generate_transcripts_from_variants.transOff = 0
        for tId, varSeq, varComb in _generate_combinations(tId, vs, list(tSeq), {}, 0):
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
                pos = v.coding[tId].tranPos + offset
                usedVs[pos] = v
                offset = _incorp.get(v.type, lambda a, b, c, d, e: e)(seq, v, tId, pos, offset)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs.copy()

                if flag:
                    for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset):
                        yield s

                # update the transcript variant id
                generate_transcripts_from_tumor_variants.transOff += 1
                pos = v.coding[tId].tranPos + offset
                usedVs[pos] = v

                offset = _incorp.get(v.type, lambda a, b, c, d, e: e)(seq, v, tId, pos, offset)

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
        if strand == REVERS:
            for flag, v in transToVar[tId]:
                v.ref = v.ref[::-1].translate(COMPLEMENT)
                v.obs = v.obs[::-1].translate(COMPLEMENT)

        vs.sort(key=lambda (isTumor, v): v.genomePos-1
                if v.type in [VariationType.FSINS, VariationType.INS]
                else v.genomePos, reverse=True)
        if not _check_for_problematic_variants(map(lambda x: x[1],vs)):
            warnings.warn("Intersecting variants found for Transcript %s"%tId)
            continue

        generate_transcripts_from_tumor_variants.transOff = 0
        for tId, varSeq, varComb in _generate_combinations(tId, vs, list(tSeq), {}, 0):
            yield Transcript("".join(varSeq), geneid, tId, _vars=varComb)


################################################################################
#        T R A N S C R I P T    = = >    P R O T E I N
################################################################################
def generate_proteins_from_transcripts(transcripts, table='Standard', stop_symbol='*', to_stop=False, cds=False):
        """
        Enables the translation from a transcript to a protein instance

        :param: list(Transcript)/Transcript - a list of transcripts to translate to proteins
        :returns: (Protein) -- the protein that corresponds to the transcript
        """

        if isinstance(transcripts, Transcript):
            transcripts = [transcripts]
        else:
            if any(not isinstance(t, Transcript) for t in transcripts):
                raise ValueError("Specified input is not of type Transcript")

        for t in transcripts:
            # translate to a protein sequence
            #if len(str(self)) % 3 != 0:
            #    raise ValueError('ERROR while translating: lenght of transcript %s is no multiple of 3, the transcript is:\n %s' % (self.transcript_id, self))

            #TODO warn if intrasequence stops - biopython warns if  % 3 != 0
            prot_seq = str(t.translate(table=table, stop_symbol=stop_symbol, to_stop=to_stop, cds=cds))

            # only transfer the non-synonymous variants to the protein as an
            # ordered dict, also translate into protein positions
            new_vars = dict()
            for pos, var in t.vars.iteritems():
                if not var.isSynonymous:
                    # prot_pos = int(math.ceil(pos/3.0)) #if pos is 0 based this is wrong
                    prot_pos = pos // 3
                    # new_vars.setdefault(pos, []).append(var) #this is obviously wrong
                    new_vars.setdefault(prot_pos, []).append(var)

            gene_id = t.gene_id
            yield Protein(prot_seq, gene_id, t.transcript_id, t, new_vars)


################################################################################
#        P R O T E I N    = = >    P E P T I D E
################################################################################

def generate_peptides_from_protein(proteins, window_size, peptides=None):
    """
    Creates all peptides for a given window size, from a given protein. The
    result is a generator.

    :param Protein protein: (list of) protein(s) from which a list of unique
                            peptides should be generated
    :param int window_size: size of peptide fragments
    :param list(Peptide) peptides: a list of peptides to update during peptide generation
                                (usa case: Adding and updating Peptides of newly generated Proteins)
    """

    def gen_peptide_info(protein):
        # Generate peptide sequences and find the variants within each
        res = []

        seq = str(protein)
        for i in xrange(len(protein)+1-window_size):
            # generate peptide fragment
            end = i+window_size
            pep_seq = seq[i:end]

            res.append((pep_seq, i))
        return res

    if isinstance(peptides, Peptide):
        peptides = [peptides]

    if peptides and any(not isinstance(p, Peptide) for p in peptides):
        raise ValueError("Specified list of Peptides contain non peptide objects")

    final_peptides = {} if peptides is None else {str(p):p for p in peptides}

    if isinstance(proteins, Protein):
        proteins = [proteins]

    if any(not isinstance(p, Protein) for p in proteins):
        raise ValueError("Input does contain non protein objects.")

    for prot in proteins:
        # generate all peptide sequences per protein:
        for (seq, pos) in gen_peptide_info(prot):

            t_id = prot.transcript_id
            if seq not in final_peptides:
                final_peptides[seq] = Peptide(seq)

            final_peptides[seq].proteins[t_id] = prot
            final_peptides[seq].proteinPos[t_id].append(pos)

    return final_peptides.values()
