# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.Generator
   :synopsis: Contains functions to transform variants to transcripts and 
              proteins to peptides.

   :Note: All internal indices start at 0!

.. moduleauthor:: schubert, walzer

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
    Incorporates a snp into the given :class:`~Fred2.Core.Transcript.Transcript` sequence (side effect!).

    :param list(char) seq: :class:`~Fred2.Core.Transcript.Transcript` sequence as a list
    :param  var: The snp variant to incorporate
    :type var: :class:`~Fred2.Core.Variant.Variant`
    :param str transId: The transcript ID of seq
    :param int pos: The position of the variant
    :param int offset: The offset which has to be added onto the :class:`~Fred2.Core.Transcript.Transcript` position of
                       the :class:`~Fred2.Core.Variant.Variant`
    :param bool isReverse: Defines whether current transcript is reverse oriented
    :return: The the modified offset
    :rtype: int
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
    Incorporates an insertion into the given :class:`~Fred2.Core.Transcript.Transcript` sequence (side effect!).

    :param list(char) seq: :class:`~Fred2.Core.Transcript.Transcript` sequence as a list
    :param var: The snp variant to incorporate
    :type var: :class:`~Fred2.Core.Variant.Variant`
    :param str transId: The :class:`~Fred2.Core.Transcript.Transcript` ID of seq
    :param int offset: The offset which has to be added onto the :class:`~Fred2.Core.Transcript.Transcript` position of
                       the :class:`~Fred2.Core.Variant.Variant`
    :param bool isReverse: Whether transcript is reverse oriented
    :return: The modified offset
    :rtype: int
    """
    if var.type not in [VariationType.INS, VariationType.FSINS]:
        raise TypeError("%s is not a insertion"%str(var))

    obs = var.obs[::-1].translate(COMPLEMENT) if isReverse else var.obs

    seq[pos:pos] = obs
    return offset + len(obs)


def _incorp_deletion(seq, var, transId, pos, offset, isReverse=False):
    """
    Incorporates a deletion into the given transcript sequence (side effect!).

    :param list(char) seq: :class:`~Fred2.Core.Transcript.Transcript` sequence as a list
    :param var: The snp variant to incorporate
    :type var: :class:`~Fred2.Core.Variant.Variant`
    :param str transId: The :class:`~Fred2.Core.Transcript.Transcript` ID of seq
    :param int offset: The offset which has to be added onto the :class:`~Fred2.Core.Transcript.Transcript` position of
                       the :class:`~Fred2.Core.Variant.Variant`
    :param bool isReverse: Whether transcript is reverse oriented
    :return: The modified offset
    :rtype: int
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
    Tests for problematic variants, e.g. variants that coincide.

    :param vars: Initial list of variants
    :type vars:  list(:class:`~Fred2.Core.Variant.Variant`)
    :return: True if no intersecting variants were found
    :rtype: bool
    :invariant: list(:class:`~Fred2.Core.Variant.Variant`) vars: List is sorted based on genome position in descending
                                                                 order
    """
    def get_range(var):
        current_range = [0, 0]
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
def generate_peptides_from_variants(vars, length, dbadapter, id_type, peptides=None,
                                    table='Standard', stop_symbol='*', to_stop=True, cds=False, db="hsapiens_gene_ensembl"):
    """
    Generates :class:`~Fred2.Core.Peptide.Peptide` from :class:`~Fred2.Core.Variant.Variant` and avoids the
    construction of all possible combinations of heterozygous variants by considering only those within the peptide
    sequence window. This reduces the number of combinations from 2^m with m = #Heterozygous Variants to 2^k with
    k<<m and k = #Heterozygous Variants within peptide window (and all frame-shift mutations that occurred prior to
    the current peptide window).

    The result is a generator.

    :param vars: A list of variant objects to construct peptides from
    :type vars: list(:class:`~Fred2.Core.Variant.Variant`)
    :param int length: The length of the peptides to construct
    :param dbadapter: A :class:`~Fred2.IO.ADBAdapter.ADBAdapter` to extract relevant transcript information
    :type dbadapter: :class:`~Fred2.IO.ADBAdapter.ADBAdapter`
    :param id_type: The type of the transcript IDs used in annotation of variants (e.g. REFSEQ, ENSAMBLE)
    :type id_type: :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`
    :param peptides: A list of pre existing peptides that should be updated
    :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`)
    :param str table: Which codon table to use? This can be either a name (string), an NCBI identifier (integer), or
                         a CodonTable object (useful for non-standard genetic codes). Defaults to the 'Standard' table
    :param str stop_symbol: Single character string, what to use for any terminators, defaults to the asterisk, '*'
    :param bool to_stop: Boolean, defaults to False meaning do a full translation continuing on past any stop codons
                        (translated as the specified stop_symbol). If True, translation is terminated at the first in
                        frame stop codon (and the stop_symbol is not appended to the returned protein sequence)
    :param bool cds: cds - Boolean, indicates this is a complete CDS. If True, this checks the sequence starts with a
                     valid alternative start codon (which will be translated as methionine, M),
                     that the sequence length is a multiple of three, and that there is a single in frame stop codon at
                     the end (this will be excluded from the protein sequence, regardless of the to_stop option). If
                     these tests fail, an exception is raised
    :return: A list of unique (polymorphic) peptides
    :rtype: Generator(:class:`~Fred2.Core.Peptide.Peptide`)
    :raises ValueError: If incorrect table argument is pasted
    :raises TranslationError: If sequence is not multiple of three, or first codon is not a start codon, or last codon
                              is not a stop codon, or an extra stop codon was found in frame, or codon is non-valid
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
                offset = _incorp.get(v.type, lambda a, b, c, d, e, f: e)(seq, v, tId, pos, offset, isReverse)
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
                offset = _incorp.get(v.type, lambda a, b, c, d, e, f: e)(seq, v, tId, pos, offset, isReverse)

                for s in _generate_combinations(tId, vs, seq, usedVs, offset, isReverse):
                    yield s
        else:
            yield tId + ":FRED2_%i"%transOff, seq, usedVs

    if not isinstance(dbadapter, ADBAdapter):
        raise TypeError("The given dbadapter is not of type ADBAdapter")

    transToVar = {}
    for v in vars:
        for trans_id in v.coding.keys():
            transToVar.setdefault(trans_id, []).append(v)

    prots = []
    for tId, vs in transToVar.items():
        #print tId
        query = dbadapter.get_transcript_information(tId, type=id_type, _db=db)
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
        for start in range(0, len(tSeq) + 1 - 3 * length, 3):
            #supoptimal as it always has to traverse the combination tree for all frameshift mutations.
            end = start + 3 * length
            vars_fs_hom = [x for x in vs if (x.isHomozygous
                                    or x.type in [VariationType.FSINS, VariationType.FSDEL])
                                    and x.coding[tId].tranPos < start]
            vars_in_window = [x for x in vs if start <= x.coding[tId].tranPos < end]

            if vars_in_window:
                vars = vars_fs_hom+vars_in_window
                for ttId, varSeq, varComb in _generate_combinations(tId, vars, list(tSeq), {}, 0, strand == REVERS):
                    prots = chain(prots, generate_proteins_from_transcripts(Transcript("".join(varSeq), geneid, ttId,
                                                                                       vars=varComb)))
    return [ p for p in generate_peptides_from_proteins(prots, length, peptides=peptides)
             if any(p.get_variants_by_protein(prot) for prot in p.proteins.keys())]

################################################################################
#        V A R I A N T S     = = >    T R A N S C R I P T S
################################################################################


def generate_transcripts_from_variants(vars, dbadapter, id_type, db="hsapiens_gene_ensembl"):
    """
    Generates all possible transcript :class:`~Fred2.Core.Transcript.Transcript` based on the given
    :class:`~Fred2.Core.Variant.Variant`.

    The result is a generator.

    :param vars: A list of variants for which transcripts should be build
    :type vars: list(:class:`~Fred2.Core.Variant.Variant`)
    :param: dbadapter: a DBAdapter to fetch the transcript sequences
    :type dbadapter: class:`~Fred2.IO.ADBAdapter.ADBAdapter`
    :param id_type: The type of the transcript IDs used in annotation of variants (e.g. REFSEQ, ENSAMBLE)
    :type id_type: :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`
    :return: A generator of transcripts with all possible variations determined by the given variant list
    :rtype: Generator(:class:`~Fred2.Core.Transcript.Transcript)
    :invariant: Variants are considered to be annotated from forward strand, regardless of the transcripts real
                orientation
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
        for trans_id in v.coding.keys():
            transToVar.setdefault(trans_id, []).append(v)

    for tId, vs in transToVar.items():
        query = dbadapter.get_transcript_information(tId, type=id_type, _db=db)
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
            yield Transcript("".join(varSeq), geneid, tId, vars=varComb)


################################################################################
#        T R A N S C R I P T    = = >    P R O T E I N
################################################################################
def generate_proteins_from_transcripts(transcripts, table='Standard', stop_symbol='*', to_stop=True, cds=False):
        """
        Enables the translation from a :class:`~Fred2.Core.Transcript.Transcript` to a
        :class:`~Fred2.Core.Protein.Protein` instance. The result is a generator.

        The result is a generator.

        :param transcripts:  A list of or a single transcripts to translate
        :type transcripts: list(:class:`~Fred2.Core.Transcript.Transcript`) or :class:`~Fred2.Core.Transcript.Transcript`
        :param str table: Which codon table to use? This can be either a name (string), an NCBI identifier (integer),
                          or a CodonTable object (useful for non-standard genetic codes). Defaults to the 'Standard'
                          table
        :param str stop_symbol: Single character string, what to use for any terminators, defaults to the asterisk, '*'
        :param bool to_stop: Translates sequence and passes any stop codons if False (default True)(translated as the
                             specified stop_symbol). If True, translation is terminated at the first in frame stop
                             codon (and the stop_symbol is not appended to the returned protein sequence)
        :param bool cds: Boolean, indicates this is a complete CDS. If True, this checks the sequence starts with a
                         valid alternative start codon (which will be translated as methionine, M), that the sequence
                         length is a multiple of three, and that there is a single in frame stop codon at the end
                         (this will be excluded from the protein sequence, regardless of the to_stop option).
                         If these tests fail, an exception is raised
        :returns: The protein that corresponds to the transcript
        :rtype: Generator(:class:`~Fred2.Core.Protein.Protein`)
        :raises ValueError: If incorrect table argument is pasted
        :raises TranslationError: If sequence is not multiple of three, or first codon is not a start codon, or last
                                  codon ist not a stop codon, or an extra stop codon was found in frame, or codon is
                                  non-valid

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
            for pos, var in t.vars.items():
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
    Creates all :class:`~Fred2.Core.Peptide.Peptide` for a given window size, from a given
    :class:`~Fred2.Core.Protein.Protein`.

    The result is a generator.

    :param proteins: (Iterable of) protein(s) from which a list of unique peptides should be generated
    :type proteins: list(:class:`~Fred2.Core.Protein.Protein`) or :class:`~Fred2.Core.Protein.Protein`
    :param int window_size: Size of peptide fragments
    :param peptides: A list of peptides to update during peptide generation (usa case: Adding and updating Peptides of
                     newly generated Proteins)
    :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`)
    :return: A unique generator of peptides
    :rtype: Generator(:class:`~Fred2.Core.Peptide.Peptide`)
    """

    def gen_peptide_info(protein):
        # Generate peptide sequences and returns the sequence
        # #and start position within the protein
        res = []

        seq = str(protein)
        for i in range(len(protein)+1-window_size):
            # generate peptide fragment
            end = i+window_size
            pep_seq = seq[i:end]
            res.append((pep_seq, i))
        return res

    if isinstance(peptides, Peptide):
        peptides = [peptides]

    final_peptides = {}

    if peptides:
        for p in peptides:
            if not isinstance(p, Peptide):
                raise ValueError("Specified list of Peptides contain non peptide objects")
            final_peptides[str(p)] = p

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

    return iter(final_peptides.values())
