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


def _incorp_snp(seq, var, transId, pos, offset, isReverse=False):
    """
    Incorporates a snp into the given transcript sequence

    :param list(char) seq: Transcript sequence as a list
    :param Variant var: The snp variant to incorporate
    :param str transId: The transcript ID of seq
    :param int pos: The position of the variant
    :param int offset: The offset which has to be added onto the transcript
                       position of the variant
    :param bool isReverse: Defines whether current transcript is reverse oriented
    :return: int - The the modified offset
    """
    if VariationType.SNP != var.type:
        raise TypeError("%s is not a SNP"%str(var))
   
    ref = var.ref[::-1].translate(COMPLEMENT) if isReverse else var.ref
    obs = var.obs[::-1].translate(COMPLEMENT) if isReverse else var.obs

    #print transId, len(seq), var.get_transcript_position(transId)-1
    if seq[pos] != ref:
        warnings.warn("For %s bp does not match ref of assigned variant %s. Pos %i, var ref %s, seq ref %s " % (
        transId, str(var), pos, ref,
        seq[pos]))

    seq[pos] = obs

    return offset


def _incorp_insertion(seq, var, transId, pos, offset, isReverse=False):
    """
    Incorporates an insertion into the given transcript sequence

    :param list(char) seq: Transcript sequence as a list
    :param Variant var: The snp variant to incorporate
    :param str transId: The transcript ID of seq
    :param int offset: The offset which has to be added onto the transcript
                       position of the variant
    :param bool isReverse: Whether transcript is reverse oriented
    :return: int - The modified offset
    """
    if var.type not in [VariationType.INS, VariationType.FSINS]:
        raise TypeError("%s is not a insertion"%str(var))

    obs = var.obs[::-1].translate(COMPLEMENT) if isReverse else var.obs

    seq[pos:pos] = obs
    return offset + len(obs)


def _incorp_deletion(seq, var, transId, pos, offset, isReverse=False):
    """
    Incorporates a deletion into the given transcript sequence

    :param list(char) seq: Transcript sequence as a list
    :param Variant var: The snp variant to incorporate
    :param str transId: The transcript ID of seq
    :param int pos: The starting position of the deletion
    :param int offset: The offset which has to be added onto the transcript
                       position of the variant
    :param bool isReverse: Whether transcript is reverse oriented
    :return: int - The modified offset
    """
    if var.type not in [VariationType.DEL, VariationType.FSDEL]:
        raise TypeError("%s is not a deletion"%str(var))

    ref = var.ref[::-1].translate(COMPLEMENT) if isReverse else var.ref

    s = slice(pos, pos+len(ref))
    del seq[s]
    return offset - len(ref)


_incorp = {
            VariationType.DEL: _incorp_deletion,
            VariationType.FSDEL: _incorp_deletion,
            VariationType.INS: _incorp_insertion,
            VariationType.FSINS: _incorp_insertion,
            VariationType.SNP: _incorp_snp
         }

_allowed_aas = frozenset('ACDEFGHIKLMNPQRSTVWY')


def _check_for_problematic_variants(vars):
    """
    Filters problematic variants, e.g. variants that coincide.

    :param list(Variant) vars: Initial list of variants
    :return: bool - True if no intersecting variants were found
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
def generate_peptides_from_variants(vars, length, dbadapter, peptides=None):
    """
    Generates peptides from variants and avoids of construction all possible combinations of heterozygious variants
    by considering only those within the peptide sequence window. This reduces the number of combinations from
    2^m with m = #Heterozygious Variants to 2^k with k<<m and k = #Henterozygious Variants within peptide window
    (and all frame-shift mutations that occurred prior to the current peptide window)

    :param list(Variant) vars: A list of variant objects to construct peptides from
    :param int length: The length of the peptides to construct
    :param ADBAdapter dbadapter: A ADBAdapter to extract relevant transcript information
    :param list(Peptide) peptides: A list of pre existing peptides that should be updated
    :return: list(Peptide) - A list of unique (polymorphic) peptides
    """

    def _generate_combinations(tId, vs, seq, usedVs, offset, isReverse):
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
                offset = _incorp.get(v.type, lambda a, b, c, d, e: e)(seq, v, tId, pos, offset, isReverse)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset,isReverse):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs.copy()

                #generate transcript without the current variant
                for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset, isReverse):
                    yield s

                #and one transcript with current variant as we can't resolve haplotypes
                # update the transcript variant id
                generate_peptides_from_variants.transOff += 1
                pos = v.coding[tId].tranPos + offset
                usedVs[pos] = v
                offset = _incorp.get(v.type, lambda a, b, c, d, e: e)(seq, v, tId, pos, offset, isReverse)

                for s in _generate_combinations(tId, vs, seq, usedVs, offset, isReverse):
                    yield s
        else:
            yield tId + ":FRED2_%i"%transOff, seq, usedVs

    if not isinstance(dbadapter, ADBAdapter):
        raise TypeError("The given dbadapter is not of type ADBAdapter")

    transToVar = {}
    for v in vars:
        for trans_id in v.coding.iterkeys():
            transToVar.setdefault(trans_id, []).append(v)

    prots = []
    for tId, vs in transToVar.iteritems():
        #print tId
        query = dbadapter.get_transcript_information(tId)
        if query is None:
            warnings.warn("Transcript with ID %s not found in DB"%tId)
            continue

        tSeq = query[EAdapterFields.SEQ]
        geneid = query[EAdapterFields.GENE]
        strand = query[EAdapterFields.STRAND]

        vs.sort(key=lambda v: v.genomePos-1
                if v.type in [VariationType.FSINS, VariationType.INS]
                else v.genomePos, reverse=True)
        if not _check_for_problematic_variants(vs):
            warnings.warn("Intersecting variants found for Transcript %s"%tId)
            continue

        generate_peptides_from_variants.transOff = 0
        for start in xrange(0, len(tSeq) + 1 - 3 * length, 3):
            #supoptimal as it always has to traverse the combination tree for all frameshift mutations.
            end = start + 3 * length
            vars_fs_hom = filter(lambda x: (x.isHomozygous
                                    or x.type in [VariationType.FSINS, VariationType.FSDEL])
                                    and x.coding[tId].tranPos < start, vs)
            vars_in_window = filter(lambda x: start <= x.coding[tId].tranPos < end, vs)

            if vars_in_window:
                vars = vars_fs_hom+vars_in_window
                for ttId, varSeq, varComb in _generate_combinations(tId, vars, list(tSeq), {}, 0, strand == REVERS):
                    prots = chain(prots, generate_proteins_from_transcripts(Transcript("".join(varSeq), geneid, ttId,
                                                                                       _vars=varComb)))
    return generate_peptides_from_proteins(prots, length, peptides=peptides)

################################################################################
#        V A R I A N T S     = = >    T R A N S C R I P T S
################################################################################


def generate_transcripts_from_variants(vars, dbadapter):
    """
    Generates all possible transcript variations of the given variants

    :param list(Variant) vars: A list of variants for which transcripts should 
                               be build
    :param ADBAdapter dbadapter: a DBAdapter to fetch the transcript sequences
    :return: Generator(Transcripts) - A generator of transcripts with all
             possible variations determined by the given
             variant list
   :invariant: Variants are considered to be annotated from forward strand, 
            regardless of the transcripts real orientation   
    """
    def _generate_combinations(tId, vs, seq, usedVs, offset, isReverse=False):
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
                offset = _incorp.get(v.type, lambda a, b, c, d, e, f: e)(seq, v, tId, pos, offset, isReverse)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset, isReverse):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs.copy()

                for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset, isReverse):
                    yield s

                # update the transcript variant id
                generate_transcripts_from_variants.transOff += 1
                pos = v.coding[tId].tranPos + offset
                usedVs[pos] = v

                offset = _incorp.get(v.type, lambda a, b, c, d, e, f: e)(seq, v, tId, pos, offset, isReverse)

                for s in _generate_combinations(tId, vs, seq, usedVs, offset, isReverse):
                    yield s
        else:
            yield tId+":FRED2_%i"%transOff, seq, usedVs

    #1) get all transcripts and sort the variants to transcripts

    #For a transcript do:
        #A) get transcript sequences
        #B) generate all possible combinations of variants
        #C) apply variants to transcript and generate transcript object

    if not isinstance(dbadapter, ADBAdapter):
        raise TypeError("The given dbadapter is not of type ADBAdapter")

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

        vs.sort(key=lambda v: v.genomePos-1
                if v.type in [VariationType.FSINS, VariationType.INS]
                else v.genomePos, reverse=True)
        if not _check_for_problematic_variants(vs):
            warnings.warn("Intersecting variants found for Transcript %s"%tId)
            continue
        generate_transcripts_from_variants.transOff = 0
        for tId, varSeq, varComb in _generate_combinations(tId, vs, list(tSeq), {}, 0, isReverse=strand == REVERS):
            yield Transcript("".join(varSeq), geneid, tId, _vars=varComb)


def generate_transcripts_from_tumor_variants(normal, tumor, dbadapter):
    """
    generates all possible transcript variations of the given variants

    :param list(Variant) normal: A list of variants of the normal tissue
    :param list(Variant) tumor: A list of variant of the cancer tissue for which transcript should be generated
    :param ADBAdapter dbadapter: a DBAdapter to fetch the transcript sequences
    :return: Generator(Transcripts) - A generator of transcripts with all
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
        raise TypeError("The given dbadapter is not of type ADBAdapter")

    transToVar = {}
    for v in tumor:
        for trans_id in v.coding.iterkeys():
            transToVar.setdefault(trans_id, []).append((False, v))

    for v in normal:
        for trans_id in v.coding.iterkeys():
            if trans_id in transToVar:
                transToVar.setdefault(trans_id, []).append((True, v))

    for tId, vs in transToVar.iteritems():
        query = dbadapter.get_transcript_information(unknown=tId)
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

        :param: list(Transcript)/Transcript - A list of transcripts to translate to proteins
        :returns: Protein - the protein that corresponds to the transcript
        """

        if isinstance(transcripts, Transcript):
            transcripts = [transcripts]

        for t in transcripts:
            if not isinstance(t, Transcript):
                raise ValueError("An element of specified input is not of type Transcript")
            # translate to a protein sequence
            #if len(str(self)) % 3 != 0:
            #    raise ValueError('ERROR while translating: lenght of transcript %s is no multiple of 3, the transcript is:\n %s' % (self.transcript_id, self))

            #TODO warn if intrasequence stops - biopython warns if  % 3 != 0
            prot_seq = str(t.translate(table=table, stop_symbol=stop_symbol, to_stop=to_stop, cds=cds))

            new_vars = dict()
            for pos, var in t.vars.iteritems():
                if not var.isSynonymous:
                    prot_pos = pos // 3
                    new_vars.setdefault(prot_pos, []).append(var)

            gene_id = t.gene_id
            yield Protein(prot_seq, gene_id, t.transcript_id, t, new_vars)


################################################################################
#        P R O T E I N    = = >    P E P T I D E
################################################################################

def generate_peptides_from_proteins(proteins, window_size, peptides=None):
    """
    Creates all peptides for a given window size, from a given protein. The
    result is a generator.

    :param Protein protein: (Iterable of) protein(s) from which a list of unique
                            peptides should be generated
    :param int window_size: Size of peptide fragments
    :param list(Peptide) peptides: A list of peptides to update during peptide generation
                                (usa case: Adding and updating Peptides of newly generated Proteins)
    :return list(Peptide) - A unique list of peptides
    """

    def gen_peptide_info(protein):
        # Generate peptide sequences and returns the sequence
        # #and start position within the protein
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

    for prot in proteins:
        if not isinstance(prot, Protein):
            raise ValueError("Input does contain non protein objects.")
        # generate all peptide sequences per protein:
        for (seq, pos) in gen_peptide_info(prot):
            if all(a in _allowed_aas for a in seq.upper()):
                t_id = prot.transcript_id
                if seq not in final_peptides:
                    final_peptides[seq] = Peptide(seq)
                final_peptides[seq].proteins[t_id] = prot
                final_peptides[seq].proteinPos[t_id].append(pos)

    return final_peptides.values()
