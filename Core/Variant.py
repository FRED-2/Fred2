# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
import logging
from Base import MetadataLogger


class MutationSyntax():
    def __init__(self):
        self.type = None
        self.cds_mutation_syntax = None
        self.aa_mutation_syntax = None


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
        self.coding = None #dict transcript:MutationSyntax

        for meta in metadata:
            self.log_metadata(meta, metadata[meta])

        # TODO try to get gene and coding from metadata with IO.parse_annovar_annotation and ?

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