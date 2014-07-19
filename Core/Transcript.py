# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
import warnings
import logging
import itertools
from operator import itemgetter, attrgetter

from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna

from Protein import Protein, ProteinSet
from Variant import Variant
from Base import MetadataLogger, FrameshiftNode


class Transcript(MetadataLogger):

    def __init__(self, transcript_id, transcript_seq, protein_id=None, protein_seq=None):
        """

        :param transcript_id: (String) Transcript RefSeqID
        :param transcript_seq: (Bio.Seq) Transcript RefSeq sequence
        :param protein_id:
        :param protein_seq:
        """
        MetadataLogger.__init__(self)
        self.id = transcript_id
        # cast our sequence into a Bio.Seq object with the correct alphabet if needed
        self.sequence = Seq(transcript_seq, generic_rna) if not isinstance(transcript_seq, Seq) else transcript_seq

        self.variants = list()

    def __len__(self):
        return len(self.sequence)

    def add_variant(self, variant):
        self.variantset.append(variant)

    def extend_variants(self, other):
        self.variantset.extend(other)

    def integrate_variants(self):
        self.variantsets.sort(key=attrgetter('start'))
        #TODO create translation

    def translate_with_fs(self, frameshifts=None):
        # frameshifts is a dict in {pos: Variant} form. NOT VariantSet! We are translating
        # with a particular FS combination and NOT calculating possible combinations here.
        if frameshifts is None:
            frameshifts = []
        else:
            frameshifts = sorted(frameshifts)  # should be already sorted, but...

        # the number of bases gained or lost by each frameshift. Positive: gain, negative: lost
        fs_shifts = [(fpos, fpos[0] - fpos[1] + len(fsvar)) for
            fpos, fsvar in frameshifts]

        def reposition(orig_pos):
            start, stop = orig_pos
            new_start, new_stop = start, stop
            for (fs_start, fs_stop), fs_shift in fs_shifts:
                if fs_start <= start < fs_stop or fs_start < stop <= fs_stop:
                    warnings.warn('Watch out, variant inside frameshift! We\'re not ready to handle '
                        'that yet. %s, (%d-%d)' % (self.id, fs_start, fs_stop))
                if start >= fs_stop:  # frameshift happened before variant, so variant shifts
                    new_start += fs_shift
                    new_stop += fs_shift
            return new_start, new_stop

        fs_positions = []
        new_seq = Seq('', generic_nucleotide)
        original_seq = self.sequence[self.cds[0]:]
        next_start = 0
        for (fs_start, fs_stop), fs_var in frameshifts:
            new_seq += original_seq[next_start:fs_start]
            fs_positions.append(len(new_seq)/3)  # register first AA position that current FS affects
            new_seq += fs_var.sequence
            next_start = fs_stop
        else:
            new_seq += original_seq[next_start:]

        protein = Protein(new_seq.translate(), self)

        # now with the new sequence created it's time to translate non-FS variants. Since the frameshifts
        # moved their relative positions around, we have to use their updated locations.

        new_variantsets = {}
        for (start, stop), vset in {reposition(vpos): vset for vpos, vset in self.variantsets.iteritems()}.iteritems():
            cstart = start - (start % 3)  # codon start
            cstop = (stop + 2) / 3 * 3  # codon stop
            new_vset = VariantSet(vset.genomic_pos, set([]))

            # TODO: this may introduce superfluous AA-s, that is 'Q'->'QP' when a
            # ''->'P' would be enough. Need to look into it. -- 99% SOLVED.
            for v in vset:
                if v.variant_type not in ('FSI', 'FSD'):
                    aa_seq = (new_seq[cstart:start] + v.sequence + new_seq[stop:cstop]).translate()
                    translated_variant = Variant(v.genomic_pos, v.variant_type, aa_seq, 'AA', v.sample_id)
                    # TODO: should we carry over metadata? I think we really should!
                    # for now, let's just keep a simple reference to the original variant
                    translated_variant.log_metadata('origin', v)
                    new_vset.add_variant(translated_variant)
                    new_vset.log_metadata('origin', vset)  # TODO: maybe origin should be a first-class attribue not metadata?
            if len(new_vset) > 0:  # frameshift VariantSets would create empty new_vsets, disregard them
                new_variantsets[(cstart/3, cstop/3)] = new_vset

        protein.variantsets = new_variantsets
        protein._trim_after_stop()

        # now let's see which frameshifts were actually kept. As induced stop codons may have terminated
        # the translated sequence, there's a chance that later frameshifts are irrelevant.

        # <= instead of < as the stop codon (Biopython '*') is trimmed away and if a FS induces that
        # as its first affected AA position it DID play a role in what the sequence has become
        # although '*' is not part of the protein sequence itself.
        fs_positions = filter(lambda x: x<=len(protein), fs_positions)
        used_frameshifts = zip(fs_positions, (fs for _, fs in frameshifts[:len(fs_positions)]))

        assert protein.get_metadata('frameshifts') == [], ("Someone has tweaked with the 'frameshift'"
            " field of protein metadata before. May have come from inherited transcript metadata."
            " Use a different field name in your custom functions.")
        protein.log_metadata('frameshifts', used_frameshifts)
        return protein

    def translate(self, with_variants=False):

        frameshifts = []
        for vpos, vset in sorted(self.variantsets.iteritems()):
            for v in vset:
                if v.variant_type in ('FSI', 'FSD'):
                    frameshifts.append((vpos, v.get_metadata('genotype')[0], v))

        node_below = FrameshiftNode(0, None)

        for fs in frameshifts[::-1]:
            genotype = fs[1]
            new_node = FrameshiftNode(fs, node_below, genotype)
            node_below = new_node

        legal_combinations = node_below.traverse()
        fs_to_generate = []
        for combination in legal_combinations:
            fs_to_generate.append([(pos, fsvar) for pos, genotype, fsvar in (c.fs for c in combination)])

        result_proteins = []
        while fs_to_generate:
            fs_combination = fs_to_generate[0]
            fs_prot = self.translate_with_fs(fs_combination)
            used_fs = [fs for _, fs in fs_prot.get_metadata('frameshifts')[0]]
            # we delete all upcoming combinations that begin with the frameshifts that we've just
            # used for the previous guy. Reason: they will lead to the exact same protein. Since
            # the combination list (thanks to the frameshift combinations' enumeration method) is
            # sorted in a way that for any i < j < k
            # [fs_1, ..., fs_i, fs_j, ...] is always before
            # [fs_1, ..., fs_i, fs_k, ...], we can safely throw away the latter if we know that
            # a stop codon terminates translation somewhere between fs_i and fs_j.
            # Outer loop: filtering, inner loop: getting 2nd fields (frameshift Variant objects)
            # from a combination
            fs_to_generate = [comb for comb in fs_to_generate if [c for _, c in comb[:len(used_fs)]]!=used_fs]
            result_proteins.append(fs_prot)

        return result_proteins

    def filter(self, sample_ids, keep=False):
        filtered_vsets = {}  # new variantsets dict
        for vpos, vset in self.variantsets.iteritems():
            filtered_vset = vset.filter(sample_ids, keep)
            if filtered_vset:  # may return None if everything was filtered out
                filtered_vsets[vpos] = filtered_vset
        if filtered_vsets:
            filtered_transcript = Transcript(self.id)  # already in pool, so this is enough
            filtered_transcript.variantsets = filtered_vsets
            return filtered_transcript
        else:
            return None

    def variant_locations(self):
        variant_positions = sorted(self.variantsets.keys())
        if len(variant_positions) < 2:
            return
        print self.id
        print '\t', '\t'.join([str(v[0]) for v in variant_positions])
        print '\t', '\t'.join([str(v[1]) for v in variant_positions])

    def frameshifts(self):
        # stupid old debug thing, will get rid of this
        fshifts = []
        for vpos, vset in self.variantsets.iteritems():
            fsvar = [v.sequence for v in vset if v.variant_type in ('FSI', 'FSD')]
            if fsvar:
                if len(fsvar[0]) % 3 == 1:
                    fshifts.append('+')
                else:
                    fshifts.append('-')
        return ''.join(fshifts)

    def sanity_check(self):
        def get_sanity_vector():
            match, complement, mismatch = 0, 0, 0
            for vpos, vset in self.variantsets.iteritems():
                vstart, vstop = vpos
                cds_start, cds_stop = self.cds
                in_seq = self.sequence[cds_start+vstart:cds_start+vstop]
                in_ref = [vv.sequence for vv in vset if vv.variant_type == "REF"][0]

                # gotta do the str() to extract sequence because Seq() doesn't suport
                # equality check yet.
                if str(in_seq) == str(in_ref):
                    match += 1
                elif str(in_seq) == str(in_ref.reverse_complement()):
                    complement += 1
                else:
                    mismatch += 1
            # in an ideal world this would always be x, 0, 0.
            return match, complement, mismatch

        match, complement, mismatch = get_sanity_vector()

        # this triplet should now be x, 0, 0 if everything was fine.
        # Reality is we still have mismatching parts though in a few percents of variants.
        # In one specific case we figured out that it was simply due to the fact that the DNA
        # reference had T and the RNA reference had G on that position, and Annovar annotated it
        # with the DNA reference and thus didn't match because we use RNA reference.
        # Since then we moved to a Transcript database that Annovar uses but sequences are still
        # not in 1-1 relation with transcript metadata information (Annovar's fault).
        return match, complement, mismatch

    def print_vars(self, wrap=80):

        wrap = wrap / 3 * 3  # round down to a multiple of 3 so codons don't wrap over lines
        raw_seq = str(self.sequence)
        raw_seq_codonize = raw_seq[:self.cds[0]].lower()
        c = self.cds[0]
        raw_seq_codonize += ''.join((raw_seq[c+i:c+i+3].upper() if i%6 == 0 else
            raw_seq[c+i:c+i+3].lower() for i in range(0, len(raw_seq)-c, 3)))
        raw_seq = raw_seq_codonize

        lines = [' ' * len(raw_seq)]

        def insert_at(position, annotation):
            span = len(annotation)
            for i, ll in enumerate(lines):
                if ll[position:position+span].strip(' ') == '':
                    # we have room to insert our annotation (all spaces)
                    lines[i] = ll[:position] + annotation + ll[position+span:]
                    break
            else:  # none of the annotation lines were blank at the required position
                lines.append(' ' * len(raw_seq))  # ...so add new blank line
                insert_at(position, annotation)  # ...and add our annotation

        insert_at(self.cds[0], '\CDS start')
        insert_at(self.cds[1]-8, 'CDS end/')

        for (vstart, vend), vset in self.variantsets.iteritems():
            for vv in vset:
                if vv.variant_type in ('INS', 'FSI'):
                    insert_at(self.cds[0] + vstart, '\%s' % str(vv.sequence))
                elif vv.variant_type in ('DEL', 'FSD'):
                    insert_at(self.cds[0] + vstart, 'x' * (vend - vstart + 1)) # TODO
                elif vv.variant_type in ('SNV', 'SNP'):
                    insert_at(self.cds[0] + vstart, str(vv.sequence))

        lines = [raw_seq] + lines
        # pad first line with spaces so that we can wrap at codon boundary
        lines = ['   '[:-(self.cds[0] % 3)] + ll for ll in lines]
        if wrap:
            for i in range(0, len(raw_seq), wrap):
                for is_annot_line, ll in enumerate(lines):
                    curr_line = ll[i:i+wrap]
                    if curr_line.strip(' ') != '':
                        print '.' if is_annot_line else ' ', curr_line
        else:
            print '\n'.join(lines)


class TranscriptSet(dict, MetadataLogger):

    def __init__(self, elements=None, refseq_db=None):

        if elements == None:
            elements = {}  # NOT in function definition! dict is mutable.
        else:
            assert all((isinstance(t, Transcript) for t in elements.itervalues())), (
                'Non-transcript element found in TranscriptSet-initializing dictionary')
        dict.__init__(self, elements)
        MetadataLogger.__init__(self)

        #elements is a refseqID
        #TODO

        #elements is a RefSeqDB.get_transcript result
        #TODO

        #elements is variant set
        genes = set([g for g.gene in elements])
        for g in genes:
            ts = TranscriptSet(refseq_db.get_transcript(g))
            ts.addVariants(elements.filter(g))
            self.update(ts)

    def add_transcript(self, transcript):
        if transcript.id not in self:
            self[transcript.id] = transcript
        else:
            warnings.warn('Failed to add transcript to TranscriptSet as a transcript with a '
                          'similar ID was already there.')

    def merge(self, other):
        for tr_id, other_transcript in other.iteritems():
            if tr_id not in self:  # new transcript record for the new guy
                self[tr_id] = other_transcript
            else:  # need to merge variants on an already existing transcript
                own_transcript = self[tr_id]
                own_transcript.extend_variantsets(other_transcript.variantsets)

    def filter_with_sampleid(self, sample_ids):
        filtered_trset = {tr_id: transcript for tr_id, transcript in self.iteritems() if
            [1 for vset in transcript.variantsets.itervalues() if  # non-empty if at least one list below is non-empty
                [1 for v in vset if v.sample_id in sample_ids]]  # non-empty if at least one variant is in sample_ids
        }

        return TranscriptSet(filtered_trset)

    # TODO: this thing doesn't seem useful now at all. Consider redesign.
    def filter(self, sample_ids, keep=False):
        filtered_trset = {}
        # right now criteria is sample_id. We can think up a nice way of filtering for any attribute.
        for tr_id, transcript in self.iteritems():
            filtered_transcript = transcript.filter(sample_ids, keep)
            if filtered_transcript:  # may return None if no variants were left after filtering
                filtered_trset[tr_id] = filtered_transcript
        return TranscriptSet(filtered_trset)

    def translate(self):
        proteins = []
        for tt, t in self.iteritems():
            # we need itertools flattening because it's a nested list if it contains frameshifts
            # TODO but for now, let's stick to a single, non-FS version
            proteins.extend(t.translate())
            # prot_from_t = list(flatten(t.translate()))
            # proteins.extend(prot_from_t[:2])  # TODO: multiple frameshifts are not handled well yet
        return ProteinSet(proteins)
