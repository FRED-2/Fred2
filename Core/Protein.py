"""
.. module:: Protein
   :synopsis: Contains the Protein Class.
.. moduleauthor:: brachvogel

"""
__author__ = 'brachvogel'

from collections import OrderedDict
from itertools import product
import re

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import MetadataLogger




class Protein(MetadataLogger, Seq):

    def __init__(self, _seq, _gene_id, _orig_transcript=None, _vars=None):
        """
        :param _seq: Bio.Seq object consiting of an IUPACProtein alphabet
        :type _sep: Seq.
        """
        
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq, IUPAC.IUPACProtein)
        if _vars is None:
            self.vars = {}
        else:
            self.vars = _vars  # {position: variantset}
        self.orig_transcript = _orig_transcript
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
        #return str(self.sequence)
        lines = []
        for vpos, vset in self.variantsets.iteritems():
            lines.append('%s-%s: '%vpos +', '.join([('%s %s %s' % (v.variant_type, v.sample_id, v.sequence)) for v in vset]))
        return self.origin.id + '\n\t' + '\n\t'.join(lines)


def generate_peptides_from_protein(proteins, window_size):
    """
    Creates all peptides for a given window size, from a given protein. The
    result is a generator.
    :param protein:
    :type protein: Fred2.Core.Protein.Protein.
    """

    def vars_in_window(protein, start, end):
        """
        Find all variants that belong (according to their position) to a 
        peptide. If the protein did not have any variants returns None
        :param start: starting position of a specific window
        :type start: int.
        :param end: end position (inclusive) of a specific window
        :type end: int.
        """
        if protein.vars is not None:
            res = OrderedDict()

            for var in protein.vars:
                t_id = protein.orig_transcript.transcript_id # transcript-var-id
                pos = var.get_protein_position(t_id)
                if start <= pos and pos <= end:
                    res[pos] = var
            return res
        else:
            return None

    def gen_peptide_info(seq):
        """
        Generate a single peptide
        """
        res = []
        for i in range(len(seq)+1-window_size):
            pep_seq = seq[i:i+window_size]
             # get the variants within that peptide:
            pep_var = vars_in_window(seq, i, i+window_size)
            res.append(pep_seq, pep_var)
        return res



    final_peptides = {} # sequence : peptide-instance

    for prot in proteins:
        # generate all peptide sequences per protein:
        for (pep, _vars) in gen_peptide_info(str(prot)):

            if pep in final_peptides:
                t_id = prot.orig_transcript.transcript_id
                final_peptides[pep].proteins[t_id] = prot
                final_peptides[pep].vars[t_id] = _vars
                final_peptides[pep].transcripts[t_id] = prot.orig_transcript
            else:
                final_peptides[pep] = Peptide(pep)
    return final_peptides.values()




    def peptides_in_window(self, w_start, length, is_recursion=0):

        def get_vars(w_start, w_end):
            vars = []
            for (start, end), mut in sorted(self.variantsets.iteritems()):
                if w_start <= start < w_end or w_start < end < w_end:
                    vars.append(((start, end), mut))
            return vars

        w_end = w_start + length

        vars = get_vars(w_start, w_end)
        combinations = product(*[[x for x in product([pos], vset)] for pos, vset in vars])
        my_peptides = []
        my_variants = []
        # TODO: store peptide position within protein.

        for comb in combinations:

            absolute_end = max(w_end, comb[-1][0][1] if comb else 0)

            # If an insertion is found at the very beginning of the window, we'll have to push the
            # window through it, as the next window wouldn't contain the transcript at all (it would
            # occur before the window and wouldn't show up among the window's variants)
            # and we'd lose peptides which partially contain the insertion at their beginning.
            # So we store this potential insertion's length to know how many steps we'll need
            # when pushing the window through it.
            ins_length_at_start = 0
            if not is_recursion and comb:
                (v_start, v_end), v = comb[0]
                if v_start == w_start:  # insertion at 1st position of window
                    v_len = len(v) + v_start - v_end
                    if v_len > 0:
                        ins_length_at_start = v_len

            curr_pep = ''
            curr_vars = []
            curr_dels = []
            curr_fs = []

            frameshifts = self.get_metadata('frameshifts', True)
            if frameshifts:
                for fs_pos, fs_var in frameshifts:
                    if fs_pos < w_start + length:
                        curr_fs.append(fs_var)

            prev_stop = w_start
            for (v_start, v_end), var in comb:
                curr_pep += self.sequence[prev_stop:v_start]
                curr_vars.extend([None]*(v_start-prev_stop))
                # we store where the variant takes place on the PEPTIDE
                curr_pep += var.sequence
                curr_vars.extend([var]*len(var))
                if len(var) < (v_end - v_start):  # or var.variant_type=='DEL'
                    curr_dels.append((len(curr_pep), var))
                prev_stop = v_end
            curr_pep += self.sequence[prev_stop:w_end]
            curr_vars.extend([None]*(w_end-prev_stop))

            if len(curr_pep) >= length:
                if not ins_length_at_start:
                    my_peptides.append(curr_pep[:length])
                    used_variants = set(curr_vars[:length]).union(
                        set(deletion for del_pos, deletion in curr_dels if 0 < del_pos < length)).union(  # todo maybe <=
                        set(curr_fs))
                    my_variants.append(used_variants.difference([None]))
                else:  # we have an insertion at the beginning, so we have to step it through
                    for ii in range(ins_length_at_start):
                        prelim_pep = curr_pep[ii:ii+length]
                        used_variants = set(curr_vars[ii:ii+length]).union(
                            set(deletion for del_pos, deletion in curr_dels if ii < del_pos < ii+length)).union(  # todo maybe <=
                            set(curr_fs))
                        if len(prelim_pep) == length:
                            my_peptides.append(prelim_pep)
                            my_variants.append(used_variants.difference([None]))
                        else:  # we're shorter because of a deletion
                            len_missing = length - len(prelim_pep)
                            overshoot = self.peptides_in_window(absolute_end, len_missing, 1)  # set recursion flag so it doesn't do the tricky starting insertion thing
                            for oo_seq, oo_vars in overshoot:
                                my_peptides.append((prelim_pep + oo_seq)[:length])
                                used_variants = used_variants.union(oo_vars)  # todo maybe <=
                                my_variants.append(used_variants.difference([None]))

            else:  # a deletion or reaching the end made this peptide shorter. Recursively generate remaining part
                if w_end < len(self.sequence):  # we can only save this if we're not at the end
                    len_missing = length - len(curr_pep)
                    overshoot = self.peptides_in_window(absolute_end, len_missing, 1)  # this is just an overshoot, don't do the tricky starting insertion thing
                    used_variants = set(curr_vars).union(set(deletion for del_pos, deletion in curr_dels))  # todo maybe <=
                    for oo_seq, oo_vars in overshoot:
                        my_peptides.append((curr_pep + oo_seq)[:length])
                        my_variants.append(used_variants.union(oo_vars).difference([None]))

        return list(set(zip(my_peptides, map(frozenset, my_variants))))  # frozenset for hashability at a later duplicate-filtering step

    def all_peptides(self, peptide_length):
        allpeps = []
        # TODO: so I hacked up an ugly way to store window position relative to transcript but it's
        # ugly as hell and I won't be able to understand it in a few weeks. i is this position FTR.
        # Actually, rewriting the duplicate checking from list(set(allpeps)) is enough as it retains order
        # even if we don't register window locations. (as PeptideSet is an OrderedDict so it doesn't mix
        # order either!)
        for i in range(len(self)):
            onewindow = self.peptides_in_window(i, peptide_length)
            onewindow = [(i, x) for x in onewindow if '*' not in x[0]]
            allpeps.extend(onewindow)

        newallpeps = []
        seen = set()
        for (i, x) in allpeps:
            if x not in seen:
                newallpeps.append((i, x))
                seen.add(x)

        #allpeps = list(set(allpeps))  # filters stuff like QQQQQQQQQ that comes from many different windows

        pepset = PeptideSet()

        for i, (pepseq, pepvars) in newallpeps:
            peptide = Peptide(pepseq, self)
            peptide.log_metadata('variants', pepvars)
            peptide.log_metadata('rel_window', i)
            pepset.add_peptide(peptide)

        return pepset

    # SOON TO BE DEPRECATED! Seeing how fast generating all peptides (regardless of having variants
    # or not) with the cleaner method above, we can get rid of this and emulate this with brute force
    # peptide generation & filtering.
    def create_peptides(self, around_samples, include_samples, peptide_length):

        result_peptides = []
        relevant_samples = set(around_samples).union(set(include_samples))

        def find_variants_in(bstart, bend, astart, aend):
            # b: piece before variant
            # a: piece after variant
            before_vsets = {(p0, p1): vset.filter_sample_ids(relevant_samples) for (p0, p1), vset in
                self.variantsets.iteritems() if bstart < p1 <= bend}
            after_vsets  = {(p0, p1): vset.filter_sample_ids(relevant_samples) for (p0, p1), vset in
                self.variantsets.iteritems() if astart <= p0 < aend}

            # remove empty variantsets
            before_vsets = {pos: vset for pos, vset in before_vsets.iteritems() if vset.variants}
            after_vsets  = {pos: vset for pos, vset in  after_vsets.iteritems() if vset.variants}
            return before_vsets, after_vsets

        def apply_variants(pstart, pend, vsets, side="left"):
            possible_sequences = []
            # Very concise and powerful one-liner for generating all possible genotype combinations
            # of any variantset list (positionally keyed dict, in fact).
            # We only need to use it for small segments though, where most often no more than
            # one variant is present, so real combinations won't even occur.
            all_combinations = product(*[[x for x in product([pos], vset)] for pos, vset in vsets.iteritems()])

            for comb in all_combinations:
                segment_sequence = []
                last_end = pstart
                for (p0, p1), variant in comb:
                    # print "current chunk: ", pstart, pend
                    # print "applying variant now: ", p0, p1, variant.sequence, variant.sample_id
                    segment_sequence.append(str(self.sequence)[last_end:p0])
                    segment_sequence.append(str(variant.sequence))
                    # print "segment as of now: ", segment_sequence
                    last_end = p1
                segment_sequence.append(str(self.sequence)[last_end:pend])

                if side == "left":
                    segment_sequence = ''.join(segment_sequence)[-(pend-pstart):]
                else:
                    assert side == "right", "apply_variants() failed, side can only be left or right"
                    segment_sequence = ''.join(segment_sequence)[:(pend-pstart)]

                if len(segment_sequence) < (pend - pstart):
                    # a deletion made the segment shorter, so we fill the missing part
                    # TODO: check for variants in the missing few amino-acids
                    missing = pend - pstart - len(segment_sequence)
                    if side == "left":
                        segment_sequence = str(self.sequence)[pstart-missing:pstart] + segment_sequence
                    else:
                        segment_sequence += str(self.sequence)[pend:pend+missing]
                possible_sequences.append(segment_sequence)
            if not all_combinations:  # if there are no variants, just return the piece w/o variants
                possible_sequences.append(str(self.sequence)[pstart:pend])  # TODO: maybe it shouldn't do it?
            return possible_sequences

        anchor_variants = []
        for pos, vset in sorted(self.variantsets.iteritems()):
            for vv in vset:
                if vv.sample_id in around_samples and vv.variant_type != 'REF':
                    anchor_variants.append((pos, vv))
        if anchor_variants:
            print self.origin.id, 'generating peptides around:', anchor_variants

        # # should be exactly the same, but shorter
        # anchor_variants = [(pos, vv) for vv in vset for pos, vset in
        #     sorted(self.variantsets.iteritems()) if vv.sample_id in include_samples]

        for (vstart, vend), v in anchor_variants:
            Lv = len(v)  # length of variant sequence
            Lp = peptide_length  # peptide window width

            # meaning of Pb: relative position of variant's first AA compared to the peptide window
            # ranges from 1-Lv (only the last AA of the variant is in window, so the variant begins
            # as many positions earlier than the window as wide the variant is, minus one)
            # ...to Lp (only the first AA of the variant is in the window. In this case its relative
            # position equals the size of the peptide window)
            for Pb in range(1-Lv, Lp)[::-1]:
                A = max(0, -Pb)  # 1st AA of variant sequence that falls in the peptide window
                B = min(Lv, Lp - Pb)  # last AA of variant sequence that falls in the peptide window + 1
                Le = max(0, Pb)  # number of AAs in window before the variant
                Lu = max(0, Lp - Pb - Lv)  # number of AAs in window after the variant

                # these are the further variants found in the regions preceding and succeeding
                # the anchor variant. Normally they are empty, but if they DO contain close
                # variants, we need to generate all possible genotypes from them.
                # So we generate their power sets and apply the variants in all possible combinations.
                before_vsets, after_vsets = find_variants_in(vstart - Le, vstart, vend, vend + Lu)
                before_seqs = apply_variants(vstart - Le, vstart, before_vsets, "left")
                after_seqs  = apply_variants(vend, vend + Lu, after_vsets, "right")

                if before_vsets or after_vsets:
                    print 'further variants found in window:', before_vsets, after_vsets

                # power set of preceding and succeeding segment (they were independent up till now)
                before_after_combined = product(before_seqs, after_seqs)
                for before, after in before_after_combined:
                    peptide_sequence = before + v.sequence[A:B] + after
                    peptide = Peptide(peptide_sequence, self)
                    if '*' not in peptide.sequence and len(peptide) == peptide_length:
                        # if peptide contains end codon or is shorter than required (runs off protein)
                        # then don't include it
                        result_peptides.append(peptide)
                    else:
                        # TODO: check if '*' was coming from a homozygous mutation and terminate further
                        # peptide generation in that case.
                        pass

                # print (self.sequence[vstart - Le:vstart] + v.sequence[A:B] +
                #     self.sequence[vend:vend + Lu], "\t", Le, A, B, Lu, Pb, "---", before_vsets,
                #     after_vsets, before_seqs, after_seqs)
            print

        return PeptideSet(result_peptides)


    def list_variants(self):
        variants = []
        for (vstart, vend), vset in self.variantsets.iteritems():
            for v in vset:
                if v.variant_type != 'REF':
                    variants.append(v.variant_type + ":" + str(self.sequence)[vstart:vend] + str(vstart) + str(v.sequence))
        return ', '.join(variants)

    def print_vars(self, wrap=80):

        raw_seq = str(self.sequence)
        lines = [' ' * len(raw_seq)]

        def insert_at(position, annotation):
            span = len(annotation)
            for i, ll in enumerate(lines):
                if ll[position:position+span].strip(' ') == '':  # we have room to insert our annotation (all spaces)
                    lines[i] = ll[:position] + annotation + ll[position+span:]
                    break
            else:
                lines.append(' ' * len(raw_seq))  # adding new blank line
                insert_at(position, annotation)  # ...and now we'll have room to add our annotation

        for (vstart, vend), vset in self.variantsets.iteritems():
            for vv in vset:
                if vv.variant_type in ('INS', 'FSI', 'DEL', 'FSD'):
                    insert_at(vstart, 'x'*(vend-vstart) + '\%s' % str(vv.sequence))
                # elif vv.variant_type in ('DEL', 'FSD'):
                #     insert_at(vstart, 'x' * len(vv))
                elif vv.variant_type in ('SNV', 'SNP'):
                    insert_at(vstart, str(vv.sequence))

        lines = [raw_seq] + lines
        if wrap:
            for i in range(0, len(raw_seq), wrap):
                for is_annot_line, ll in enumerate(lines):
                    curr_line = ll[i:i+wrap]
                    if curr_line.strip(' ') != '':
                        print '.' if is_annot_line else ' ', curr_line
                #print
        else:
            print '\n'.join(lines)

    def origin_filter(self, attrib, attrib_values):
        return ProteinSet((protein for protein in self.proteins if
            getattr(protein.origin, attrib) in attrib_values))

    def _trim_after_stop(self):
        # in order to handle stoploss mutations and frameshifts, we don't stop translating at
        # the transcript's stop codon in order to allow translation continuing over the CDS.
        # However, if the stop codon is not overshadowed by any stoploss mutation, we DO have to
        # trim off the remaining part as it can never get translated. So we search for the first
        # unmutated stop codon and delete the sequence and variants after that.
        # This is especially important for frameshifted proteins as we have no information
        # where the first stop codon should be found.

        for a in re.finditer('\*', str(self.sequence)):  # use biopython HasStopCodon object?
            stopcodon = a.start()
            for (vstart, vstop), vset in self.variantsets.iteritems():
                if vstart <= stopcodon < vstop:
                    # variant at stop codon position! TODO: check if it is really stoploss.
                    # An insertion for example can still contain a stop codon and is thus not stoploss.
                    found_stoploss = False
                    for v in vset:
                        if '*' not in v.sequence:
                            found_stoploss = True
                            print 'Stoploss variant found!', self.origin.id
                            break
                    if found_stoploss:
                        break  # breaking this loop will mean that we'll have to continue seeking for stop codons.
            else:  # if we didn't break the previous loop but exited normally it IS a stop codon.
                self.sequence = self.sequence[:stopcodon]
                for (vstart, vstop), vset in self.variantsets.items():  # NOT iteritems...
                    # ...because we can't delete while iterating through the dict.
                    # .items() creates a copy so it's safe to delete using its values.
                    if stopcodon < vstart:
                        del self.variantsets[(vstart, vstop)]
