# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'szolek', 'walzer'

import warnings
import logging
import itertools
import re
from operator import itemgetter, attrgetter

from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
from Bio.Alphabet import generic_rna, generic_protein

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
        self.pid = protein_id
        # cast our sequence into a Bio.Seq object with the correct alphabet if needed
        self.mrna = Seq(re.sub(r'\s+', '', transcript_seq), generic_rna) if not isinstance(transcript_seq, Seq) else transcript_seq
        self.protein = Seq(re.sub(r'\s+', '', protein_seq), generic_protein) if not isinstance(protein_seq, Seq) else protein_seq

        self.variants = list()

    def __len__(self):
        return len(self.protein)

    def add_variant(self, variant):
        self.variants.append(variant)

    def extend_variants(self, other):
        self.variants.extend(other)

    def to_proteins(self, with_fs=False):
        # TODO fs
        # TODO non fs in/del
        self.variants.sort(key=attrgetter('start'))
        homos = [id(x) for x in self.variants if not bool(re.match(r"heterozygous", x.zygosity, flags=re.IGNORECASE))]
        heteros = [id(x) for x in self.variants if bool(re.match(r"heterozygous", x.zygosity, flags=re.IGNORECASE))]
        combis = list()
        for i in range(len(heteros)):
            combis += itertools.combinations(heteros, i+1)
        combis = [list(t)+homos for t in combis]
        generators = [(x for x in self.variants if id(x) in c) for c in combis]
        ret = list()
        if self.protein:
            for varcomb in generators:
                nuseq = SeqRecord(Seq(self.protein, generic_protein))
                nuseq.id = 'fred2|' + self.pid
                nuseq.name = nuseq.id
                nuseq.description = self.pid + " from " + self.id + " "
                offset = 0
                for var in varcomb:
                    if var.coding[self.id] and (var.coding[self.id].type == 'SNV'
                                                or bool(re.match(r"\Ap\.[A-X]\d+[A-X]\Z", var.coding[self.id].aa_mutation_syntax))):
                        ori, rep = re.compile(r"\d+").split(var.coding[self.id].aa_mutation_syntax.strip("p."))
                        pos = int(re.findall(r"\d+", var.coding[self.id].aa_mutation_syntax)[0])
                        if nuseq[pos-1] == ori:
                            nuseq.seq = nuseq.seq[:pos-1] + rep + nuseq.seq[pos:]
                            nuseq.description += str(var)
                        else:
                            logging.warn('No position/aminoacid match, coding is off.')
                    else:
                        logging.warn('No coding for this transcript')
                ret.append(nuseq)
        else:
            logging.warn('No protein sequence.')
        return ret

    def to_peptides(self):
        pass # TODO track variant positions

    # def translate_with_fs(self, frameshifts=None):
    #     # frameshifts is a dict in {pos: Variant} form. NOT VariantSet! We are translating
    #     # with a particular FS combination and NOT calculating possible combinations here.
    #     if frameshifts is None:
    #         frameshifts = []
    #     else:
    #         frameshifts = sorted(frameshifts)  # should be already sorted, but...
    #
    #     # the number of bases gained or lost by each frameshift. Positive: gain, negative: lost
    #     fs_shifts = [(fpos, fpos[0] - fpos[1] + len(fsvar)) for
    #         fpos, fsvar in frameshifts]
    #
    #     def reposition(orig_pos):
    #         start, stop = orig_pos
    #         new_start, new_stop = start, stop
    #         for (fs_start, fs_stop), fs_shift in fs_shifts:
    #             if fs_start <= start < fs_stop or fs_start < stop <= fs_stop:
    #                 warnings.warn('Watch out, variant inside frameshift! We\'re not ready to handle '
    #                     'that yet. %s, (%d-%d)' % (self.id, fs_start, fs_stop))
    #             if start >= fs_stop:  # frameshift happened before variant, so variant shifts
    #                 new_start += fs_shift
    #                 new_stop += fs_shift
    #         return new_start, new_stop
    #
    #     fs_positions = []
    #     new_seq = Seq('', generic_nucleotide)
    #     original_seq = self.sequence[self.cds[0]:]
    #     next_start = 0
    #     for (fs_start, fs_stop), fs_var in frameshifts:
    #         new_seq += original_seq[next_start:fs_start]
    #         fs_positions.append(len(new_seq)/3)  # register first AA position that current FS affects
    #         new_seq += fs_var.sequence
    #         next_start = fs_stop
    #     else:
    #         new_seq += original_seq[next_start:]
    #
    #     protein = Protein(new_seq.translate(), self)
    #
    #     # now with the new sequence created it's time to translate non-FS variants. Since the frameshifts
    #     # moved their relative positions around, we have to use their updated locations.
    #
    #     new_variantsets = {}
    #     for (start, stop), vset in {reposition(vpos): vset for vpos, vset in self.variantsets.iteritems()}.iteritems():
    #         cstart = start - (start % 3)  # codon start
    #         cstop = (stop + 2) / 3 * 3  # codon stop
    #         new_vset = VariantSet(vset.genomic_pos, set([]))
    #
    #         # TODO: this may introduce superfluous AA-s, that is 'Q'->'QP' when a
    #         # ''->'P' would be enough. Need to look into it. -- 99% SOLVED.
    #         for v in vset:
    #             if v.variant_type not in ('FSI', 'FSD'):
    #                 aa_seq = (new_seq[cstart:start] + v.sequence + new_seq[stop:cstop]).translate()
    #                 translated_variant = Variant(v.genomic_pos, v.variant_type, aa_seq, 'AA', v.sample_id)
    #                 # TODO: should we carry over metadata? I think we really should!
    #                 # for now, let's just keep a simple reference to the original variant
    #                 translated_variant.log_metadata('origin', v)
    #                 new_vset.add_variant(translated_variant)
    #                 new_vset.log_metadata('origin', vset)  # TODO: maybe origin should be a first-class attribue not metadata?
    #         if len(new_vset) > 0:  # frameshift VariantSets would create empty new_vsets, disregard them
    #             new_variantsets[(cstart/3, cstop/3)] = new_vset
    #
    #     protein.variantsets = new_variantsets
    #     protein._trim_after_stop()
    #
    #     # now let's see which frameshifts were actually kept. As induced stop codons may have terminated
    #     # the translated sequence, there's a chance that later frameshifts are irrelevant.
    #
    #     # <= instead of < as the stop codon (Biopython '*') is trimmed away and if a FS induces that
    #     # as its first affected AA position it DID play a role in what the sequence has become
    #     # although '*' is not part of the protein sequence itself.
    #     fs_positions = filter(lambda x: x<=len(protein), fs_positions)
    #     used_frameshifts = zip(fs_positions, (fs for _, fs in frameshifts[:len(fs_positions)]))
    #
    #     assert protein.get_metadata('frameshifts') == [], ("Someone has tweaked with the 'frameshift'"
    #         " field of protein metadata before. May have come from inherited transcript metadata."
    #         " Use a different field name in your custom functions.")
    #     protein.log_metadata('frameshifts', used_frameshifts)
    #     return protein

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
