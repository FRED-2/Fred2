# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'szolek', 'walzer'

import logging
import re
from Base import MetadataLogger


def parse_annotation(annotations, details=None):
    """
    parses an annotation string (http://www.hgvs.org/mutnomen/quickref.html)
        - e.g. from ANNOVAR: (string) ANK3:NM_001204403:exon10:c.939delG:p.M313fs,ANK3:NM_020987:exon9:c.957delG:p.M319fs,ANK3:NM_001204404:exon9:c.906delG:p.M302fs,
        - e.g. cosmic: (list of dict) ITGAV	ENST00000261023	71212	c.1502C>A	p.S501Y	...
    :param annotations: the ANNOVAR annotation string
    :param details: the ANNOVAR variant details
    :return:
    """
    mut_syn = dict()
    #gss - old ANNOVAR stuff
    # HAT1:NM_003642:exon11:c.A1208T:p.Q403L,
    # NBPF10:NM_001039703:exon4:p.A179D:GCT->GAT:+,

    #gsvar - ANNOVAR uses UCSC?
    # RANBP2:NM_006267:exon14:c.A1954G:p.I652V,
    # ANAPC1:NM_022662:exon12:c.C1393T:p.Q465X,
    # NCBP1:NM_002486:exon9:c.913_914insG:p.G305fs,
    # ANK3:NM_001204403:exon10:c.939delG:p.M313fs,ANK3:NM_020987:exon9:c.957delG:p.M319fs,ANK3:NM_001204404:exon9:c.906delG:p.M302fs,

    #cosmic - http://www.hgvs.org/mutnomen/ (http://www.hgvs.org/mutnomen/quickref.html)
    # Gene Name	Accession Number    CDS Mutation Syntax	AA Mutation Syntax
    # TP53	ENST00000269305	c.?	p.E286K
    # MAP3K6	NM_004672	c.2483C>T	p.T828I	...
    # ITGAV	ENST00000261023	c.1502C>A	p.S501Y	...
    # TP53	ENST00000269305	c.847_848insT	p.R283fs*23

    #maybe use HGVS 0.3.3dev-d13418a64c51 ??
    #TODO needs indel handling
    if isinstance(annotations, str):
        annotations = [x for x in annotations.split(',') if x]  # removes trailing list entry ''
        for a in annotations:
                    t, c, p = [None]*3
                    for i in a.split(':'):
                        if i.startswith(('NM_', 'ENST')):
                                t = i
                        elif i.startswith('c.'):
                                c = i
                        elif i.startswith('p.'):
                                p = i
                    if t and c and p:
                        mut_syn[t] = MutationSyntax(details, c, p)
                    else:
                        logging.warn('Insufficient variation coding annotation data (entry is ' + str(a) + ')')

    elif isinstance(annotations, list) \
            and all(isinstance(x, dict)
                    and 'Accession Number' in x and 'CDS Mutation Syntax' in x and 'AA Mutation Syntax' in x
                    for x in annotations):
        for a in annotations:
            mut_syn[a['Accession Number']] = MutationSyntax(details, a['CDS Mutation Syntax'], a['AA Mutation Syntax'])
    else:
        logging.warn('Unusable variation annotation input')
    return mut_syn


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
        self.observed = obs
        self.sample_id = sample_id  # hg19 if not observed in sample, but kept as reference
        self.gene = None
        self.type = None
        self.zygosity = None
        self.coding = dict()  # dict transcript_id:MutationSyntax
        self.specific = True  #TODO set variants from germline all to specific=False

        for meta in metadata:
            self.log_metadata(meta, metadata[meta])

    def __repr__(self):  # TODO, just to have something for now
        return ' '.join([self.sample_id, 'c.', self.start, self.reference, '>',  self.observed])

    def __str__(self):  # TODO this is used e.g. in transcript.to_protein for the SeqRecord description - needs pos
        return ''.join(['c.', self.start, self.reference, '>',  self.observed])  #

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
        elif 'gene' in self.metadata and self.metadata['gene'][0]:
            self.gene = self.metadata['gene'][0]
        elif refseq_DB:
            print self.gene
            self.gene = refseq_DB.get_variant_gene(self.chromosome, self.start, self.stop)
            print self.gene
        else:
            logging.warning('No gene available')
        return self.gene

    def find_coding(self, refseq_DB=None):
        for meta in self.metadata:
            if meta == 'coding' or meta == 'coding_and_splicing_details':
                self.coding.update(parse_annotation(self.metadata[meta][0]))

    def find_zygosity(self):
        if self.zygosity:
            pass
        elif 'Zygosity' in self.metadata:
            if bool(re.match(r"\Ah(omo|etero)zygous\Z", self.metadata['Zygosity'][0], flags=re.IGNORECASE)):
                self.zygosity = self.metadata['Zygosity'][0]
            else:
                self.zygosity = 'homozygous'
        elif 'genotype' in self.metadata:
            if bool(re.match(r"\Ah(omo|etero)zygous\Z", self.metadata['genotype'][0], flags=re.IGNORECASE)):
                self.zygosity = self.metadata['genotype'][0]
            else:
                self.zygosity = 'homozygous'
        else:
            logging.warning('No zygosity information available')
        return self.gene
