"""
.. module:: Transcript
   :synopsis: Contains the Transcript Class.
.. moduleauthor:: brachvogel, szolek, walzer

"""
__author__ = 'brachvogel', 'szolek', 'walzer'

#import warnings
# import logging
# import itertools
# import re
#from operator import itemgetter, attrgetter

from Bio.Seq import Seq
#from Bio.SeqIO import SeqRecord
#from Bio.Alphabet import generic_rna, generic_protein

from Fred2.Core.Protein import Protein
#from Fred2.Core.Variant import Variant
from Fred2.Core.Base import MetadataLogger
from Bio.Alphabet import IUPAC

class Transcript(MetadataLogger, Seq):
    """Transcript Class
    """

    def __init__(self, _gene_id, _transcript_id, _seq, _vars=None):
        """

        :param _transcript_id: Transcript RefSeqID
        :type _transcript_id: str.
        :param _seq: Transcript RefSeq sequence
        :type _seq: str.
        :param _start: true start position of transcript, positions starting at
                       1
        :type _start: int.
        :param _stop: true stop position of transcript
        :type _stop: int.
        :param _vars: Variants belonging to the 
            transcript
        :type _vars: {int:Fred2.Core.Variant}.
        """
        MetadataLogger.__init__(self)
        Seq.__init__(self, _seq, IUPAC.IUPACUnambiguousRNA)
        self.gene_id = _gene_id
        self.transcript_id = _transcript_id
        if vars is not None:
            self.vars = _vars
        else:
            self.vars = {}


    def __getitem__(self, index):
        """
        Overrides Bio.Seq.__getitem__

        :param index: (int) position 
        :returns: Transcript -- A Transcript consisting of the single letter at
            position :attr:`index`.
        """
        item = self[index]
        return Transcript(self.transcript_id, item, self.vars)


    def translate(self, table='Standard', stop_symbol='*', to_stop=False, 
                  cds=False):
        """
        Override of Bio.Seq.translate()
        """

        # translate to a protein sequence
        prot_seq = str(self.translate())

        # only transfer the non-synonymous variants to the protein
        new_vars = [x for x in self.vars if x.isSynonymous]

        gene_id = self.gene_id
        return Protein(prot_seq, gene_id, self, new_vars)


    # def frameshifts(self):
    #     # stupid old debug thing, will get rid of this
    #     fshifts = []
    #     for vpos, vset in self.variantsets.iteritems():
    #         fsvar = [v.sequence for v in vset if v.variant_type in ('FSI', 'FSD')]
    #         if fsvar:
    #             if len(fsvar[0]) % 3 == 1:
    #                 fshifts.append('+')
    #             else:
    #                 fshifts.append('-')
    #     return ''.join(fshifts)

    # def sanity_check(self):
    #     def get_sanity_vector():
    #         match, complement, mismatch = 0, 0, 0
    #         for vpos, vset in self.variantsets.iteritems():
    #             vstart, vstop = vpos
    #             cds_start, cds_stop = self.cds
    #             in_seq = self.sequence[cds_start+vstart:cds_start+vstop]
    #             in_ref = [vv.sequence for vv in vset if vv.variant_type == "REF"][0]

    #             # gotta do the str() to extract sequence because Seq() doesn't suport
    #             # equality check yet.
    #             if str(in_seq) == str(in_ref):
    #                 match += 1
    #             elif str(in_seq) == str(in_ref.reverse_complement()):
    #                 complement += 1
    #             else:
    #                 mismatch += 1
    #         # in an ideal world this would always be x, 0, 0.
    #         return match, complement, mismatch

    #     match, complement, mismatch = get_sanity_vector()

    #     # this triplet should now be x, 0, 0 if everything was fine.
    #     # Reality is we still have mismatching parts though in a few percents of variants.
    #     # In one specific case we figured out that it was simply due to the fact that the DNA
    #     # reference had T and the RNA reference had G on that position, and Annovar annotated it
    #     # with the DNA reference and thus didn't match because we use RNA reference.
    #     # Since then we moved to a Transcript database that Annovar uses but sequences are still
    #     # not in 1-1 relation with transcript metadata information (Annovar's fault).
    #     return match, complement, mismatch

    # def print_vars(self, wrap=80):

    #     wrap = wrap / 3 * 3  # round down to a multiple of 3 so codons don't wrap over lines
    #     raw_seq = str(self.sequence)
    #     raw_seq_codonize = raw_seq[:self.cds[0]].lower()
    #     c = self.cds[0]
    #     raw_seq_codonize += ''.join((raw_seq[c+i:c+i+3].upper() if i%6 == 0 else
    #         raw_seq[c+i:c+i+3].lower() for i in range(0, len(raw_seq)-c, 3)))
    #     raw_seq = raw_seq_codonize

    #     lines = [' ' * len(raw_seq)]

    #     def insert_at(position, annotation):
    #         span = len(annotation)
    #         for i, ll in enumerate(lines):
    #             if ll[position:position+span].strip(' ') == '':
    #                 # we have room to insert our annotation (all spaces)
    #                 lines[i] = ll[:position] + annotation + ll[position+span:]
    #                 break
    #         else:  # none of the annotation lines were blank at the required position
    #             lines.append(' ' * len(raw_seq))  # ...so add new blank line
    #             insert_at(position, annotation)  # ...and add our annotation

    #     insert_at(self.cds[0], '\CDS start')
    #     insert_at(self.cds[1]-8, 'CDS end/')

    #     for (vstart, vend), vset in self.variantsets.iteritems():
    #         for vv in vset:
    #             if vv.variant_type in ('INS', 'FSI'):
    #                 insert_at(self.cds[0] + vstart, '\%s' % str(vv.sequence))
    #             elif vv.variant_type in ('DEL', 'FSD'):
    #                 insert_at(self.cds[0] + vstart, 'x' * (vend - vstart + 1)) # TODO
    #             elif vv.variant_type in ('SNV', 'SNP'):
    #                 insert_at(self.cds[0] + vstart, str(vv.sequence))

    #     lines = [raw_seq] + lines
    #     # pad first line with spaces so that we can wrap at codon boundary
    #     lines = ['   '[:-(self.cds[0] % 3)] + ll for ll in lines]
    #     if wrap:
    #         for i in range(0, len(raw_seq), wrap):
    #             for is_annot_line, ll in enumerate(lines):
    #                 curr_line = ll[i:i+wrap]
    #                 if curr_line.strip(' ') != '':
    #                     print '.' if is_annot_line else ' ', curr_line
    #     else:
    #         print '\n'.join(lines)




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
