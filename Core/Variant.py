# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide, generic_protein
from Base import MetadataLogger


class Variant(MetadataLogger):

    def __init__(self, genomic_pos, variant_type, sequence, sequence_type, sample_id='hg19'):
        MetadataLogger.__init__(self)
        self.genomic_pos = genomic_pos  # (chr, start, end)
        self.variant_type = variant_type  # SNV, INS, DEL, FSI, FSD, REF. Maybe needed: Sec, Pyl
        self.sequence_type = sequence_type  # N, AA. Maybe sequence alphabet should suffice?
        self.sample_id = sample_id  # hg19 if not observed in sample, but kept as reference

        assert sequence_type in ('N', 'AA'), 'invalid sequence type (N or AA)'
        alphabet = generic_nucleotide if sequence_type == 'N' else generic_protein
        # make sure to cast our sequence into a Bio.Seq object with the correct alphabet if needed
        self.sequence = Seq(sequence, alphabet) if not isinstance(sequence, Seq) else sequence

    def __repr__(self):  # TODO, just to have something for now
        return ' '.join([self.variant_type, str(self.sequence), self.sample_id, self.metadata['genotype'][0]])

    def __len__(self):
        return len(self.sequence)

    # equality check does NOT take sample_id into account and it's a feature.
    # Abscract filter class will make it disappear soon though.
    def _equals(self, other):
        # BioPython raises a warning for Seq() equality comparisons but plans to
        # implement the feature. For now they advise to compare them as
        # str(seq1)==str(seq2) [and maybe the alphabets too]
        if (self.genomic_pos == other.genomic_pos
            and self.variant_type == other.variant_type
            and str(self.sequence) == str(other.sequence)
            and self.sequence.alphabet == other.sequence.alphabet):
            return True
        else:
            return False


class VariantSet(MetadataLogger):

    def __init__(self, genomic_pos, variants):
        MetadataLogger.__init__(self)
        self.genomic_pos = genomic_pos
        self.variants = set(variants)
        # these genomic_pos asserts wouldn't be needed if Variant didn't have that attribute.
        assert all((v.genomic_pos == self.genomic_pos for v in variants))

    def add_variant(self, variant):
        assert variant.genomic_pos == self.genomic_pos, 'Variant cannot be added on a different genomic position than its own'
        self.variants.add(variant)

    def filter_sample_ids(self, sample_ids):
        filtered_variants = {v for v in self.variants if v.sample_id in sample_ids}
        return VariantSet(self.genomic_pos, filtered_variants)

    # keep: if True, keeps only filtered elements, if False, throws them away
    # TODO: ugly, bad and currently unused. Will sort out in a filter redesign.
    def filter(self, sample_ids, keep=False):
        from_samples = {v for v in self.variants if v.sample_id in sample_ids}

        # all elements that are equal to anything in from_samples except for their sample_id
        matches = {v for v in self.variants for m in from_samples if v._equals(m)}
        filtered = matches if keep else self.variants.difference(matches)

        # if the filtered set would be empty, don't bother creating an empty VariantSet.
        # return None in these cases is expected behavior for higher level filter calls.
        return VariantSet(self.genomic_pos, filtered) if filtered else None

    def merge(self, otherset):  # in-place merge.
        assert self.genomic_pos == otherset.genomic_pos, 'cannot merge VariantSets on different genomic positions'

        # if both varsets contain an hg19 sample_id, they mustn't be duplicated.
        # maybe check for uniqueness of variants and drop duplicates if neccessary?
        # ...hell, even write a custom set that checks for deep equality with __hash__ interface.
        # In a normal use case, the only possible collisions could be the hg19s though.
        self.variants.update({v for v in otherset.variants if v.sample_id != 'hg19'})
        # or union and then return VariantSet(self.genomic_pos, merged)

    def __repr__(self):  # just for now. TODO
        return ', '.join([('%s %s %s' % (v.sample_id, v.variant_type, v.sequence)) for v in self.variants])

    # making VariantSet objects directly iterable, no need to iterate over vset.variants.
    # Consider subclassing from set
    def __iter__(self):
        return self.variants.__iter__()

    def __len__(self):
        return len(self.variants)
