# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'szolek', 'walzer'

import logging
import re
from Base import MetadataLogger
import IO  # TODO fix imports! (see IO)


class MutationSyntax():
    def __init__(self, type, cds, aas):
        self.type = type  # SNV, INS, DEL, FSI, FSD, REF. Maybe needed: Sec, Pyl
        self.cds_mutation_syntax = cds  #c. ...
        self.aa_mutation_syntax = aas  #p. ...


class Variant(MetadataLogger):
    def __init__(self, chrom_loc, loc_start, loc_stop, ref, obs, sample_id='hg19', metadata={}):
        MetadataLogger.__init__(self)
        self.chromosome = chrom_loc
        self.start = loc_start
        self.stop = loc_stop
        self.reference = ref
        self.observed = loc_stop
        self.sample_id = sample_id  # hg19 if not observed in sample, but kept as reference
        self.gene = None
        self.type = None
        self.zygosity = None
        self.coding = dict()  # dict transcript_id:MutationSyntax

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

    #a list of variants vl could be sorted like that vl.sort(key=lambda x: x.gene, reverse=True)
    #or sorted(vl, key=attrgetter('gene')) or for multiple criteria see http://stackoverflow.com/a/1516449

    def _equals(self, other):
        # BioPython raises a warning for Seq() equality comparisons but plans to
        # implement the feature. For now they advise to compare them as
        # str(seq1)==str(seq2) [and maybe the alphabets too]
        if (self.chromosome == other.chromosome
            and self.start == other.start
            and self.stop == other.stop
            and self.observed == other.observed
            and self.reference == other.reference
            and self.sample_id == other.sample_id):
            return True
        else:
            return False

    def get_transcript_ids(self):
        return self.coding.keys()

    def find_gene(self, refseq_DB=None):
        if self.gene:
            pass
        elif 'gene' in self.metadata:
            self.gene = self.metadata['gene'][0]
        elif refseq_DB:
            self.gene = refseq_DB.get_variant_gene(self.chromosome, self.start, self.stop)[0]
        else:
            logging.info('No gene available')
        return self.gene

    def find_coding(self, refseq_DB=None):
        for meta in self.metadata:
            if meta == 'coding' or meta == 'coding_and_splicing_details':
                self.coding.update(IO.IO.parse_annotation(self.metadata[meta][0]))

    def find_zygosity(self):
        if self.zygosity:
            pass
        elif 'Zygosity' in self.metadata \
                and bool(re.match(r"\Ah(omo|etero)zygous\Z", self.metadata['Zygosity'][0], flags=re.IGNORECASE)):
            self.zygosity = self.metadata['Zygosity'][0]
        elif 'genotype' in self.metadata \
                and bool(re.match(r"\Ah(omo|etero)zygous\Z", self.metadata['genotype'][0], flags=re.IGNORECASE)):
            self.zygosity = self.metadata['genotype'][0]
        else:
            logging.info('No zygosity information available')
        return self.gene
