# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
import logging
from Base import MetadataLogger


class Variant(MetadataLogger):

    def __init__(self, chrom_loc, loc_start, loc_stop, variant_type, ref, obs, sample_id='hg19', metadata={}):
        MetadataLogger.__init__(self)
        self.chromosome = chrom_loc
        self.start = loc_start
        self.stop = loc_stop
        self.reference = ref
        self.observed = loc_stop
        self.type = variant_type  # SNV, INS, DEL, FSI, FSD, REF. Maybe needed: Sec, Pyl
        self.sample_id = sample_id  # hg19 if not observed in sample, but kept as reference
        self.gene = None

        for meta in metadata:
            self.log_metadata(meta, metadata[meta])

    def __repr__(self):  # TODO, just to have something for now
        return ' '.join([self.type, self.sample_id, self.metadata['genotype'][0]])

    def __str__(self):  # TODO, just to have something for now
        return ' '.join([self.type, self.sample_id, self.metadata['genotype'][0]])

    def __len__(self):
        return len(self.observed)

    def __eq__(self, other):
        return self.gene == other.gene

    def __lt__(self, other):
        return self.gene < other.gene

    #vl.sort(key=lambda x: x.gene, reverse=True) a list of variants vl could be sorted like that or
    #sorted(vl, key=attrgetter('age')) or for multiple criteria http://stackoverflow.com/a/1516449

    # equality check does NOT take sample_id into account and it's a feature.
    # Abscract filter class will make it disappear soon though.
    def _equals(self, other):
        # BioPython raises a warning for Seq() equality comparisons but plans to
        # implement the feature. For now they advise to compare them as
        # str(seq1)==str(seq2) [and maybe the alphabets too]
        if (self.chromosome == other.chromosome
            and self.start == other.location_start
            and self.stop == other.location_stop
            and self.observed == other.observed
            and self.type == other.variant_type
            and self.sample_id == other.sample_id):
            return True
        else:
            return False

    #move to refseqdb class and toolbox!
    def find_gene(self, refseq_DB=None):
        if 'gene' in self.metadata:
            self.gene = self.metadata['gene'][0]
        elif refseq_DB:
            self.gene = refseq_DB.get_variant_gene(self.chromosome, self.start, self.stop)[0]
        else:
            logging.info('No gene available')

# class VariantSet(MetadataLogger):
#
#     def __init__(self, genomic_pos, variants):
#         MetadataLogger.__init__(self)
#         self.genomic_pos = genomic_pos
#         self.variants = set(variants)
#         # these genomic_pos asserts wouldn't be needed if Variant didn't have that attribute.
#         assert all((v.genomic_pos == self.genomic_pos for v in variants))
#
#     def add_variant(self, variant):
#         #this makes no sense
#         assert variant.genomic_pos == self.genomic_pos, 'Variant cannot be added on a different genomic position than its own'
#         self.variants.add(variant)
#
#     def filter_sample_ids(self, sample_ids):
#         filtered_variants = {v for v in self.variants if v.sample_id in sample_ids}
#         return VariantSet(self.genomic_pos, filtered_variants)
#
#     # keep: if True, keeps only filtered elements, if False, throws them away
#     # TODO: ugly, bad and currently unused. Will sort out in a filter redesign.
#     def filter(self, sample_ids, keep=False):
#         from_samples = {v for v in self.variants if v.sample_id in sample_ids}
#
#         # all elements that are equal to anything in from_samples except for their sample_id
#         matches = {v for v in self.variants for m in from_samples if v._equals(m)}
#         filtered = matches if keep else self.variants.difference(matches)
#
#         # if the filtered set would be empty, don't bother creating an empty VariantSet.
#         # return None in these cases is expected behavior for higher level filter calls.
#         return VariantSet(self.genomic_pos, filtered) if filtered else None
#
#     def merge(self, otherset):  # in-place merge.
#         assert self.genomic_pos == otherset.genomic_pos, 'cannot merge VariantSets on different genomic positions'
#
#         # if both varsets contain an hg19 sample_id, they mustn't be duplicated.
#         # maybe check for uniqueness of variants and drop duplicates if neccessary?
#         # ...hell, even write a custom set that checks for deep equality with __hash__ interface.
#         # In a normal use case, the only possible collisions could be the hg19s though.
#         self.variants.update({v for v in otherset.variants if v.sample_id != 'hg19'})
#         # or union and then return VariantSet(self.genomic_pos, merged)
#
#     def __repr__(self):  # just for now. TODO
#         return ', '.join([('%s %s %s' % (v.sample_id, v.variant_type, v.sequence)) for v in self.variants])
#
#     # making VariantSet objects directly iterable, no need to iterate over vset.variants.
#     # Consider subclassing from set
#     def __iter__(self):
#         return self.variants.__iter__()
#
#     def __len__(self):
#         return len(self.variants)
