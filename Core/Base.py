# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'szolek', 'walzer'

import re
from collections import defaultdict
from operator import attrgetter

from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
from Bio.Alphabet import IUPAC


class MetadataLogger(object):
    # This class provides a simple interface for assigning additional metadata to any object in our
    # data model. Examples: storing ANNOVAR columns like depth, base count, dbSNP id, quality
    # information for variants, prediction scores for peptides etc.
    # Normally used by custom toolbox functions and importers.

    def __init__(self):
        self.metadata = defaultdict(list)

    def log_metadata(self, label, value):
        self.metadata[label].append(value)

    def get_metadata(self, label, only_first=False):
        # although defaultdict *would* return [] if it didn't find label in self.metadata, it would
        # come with the side effect of adding label as a key to the defaultdict, so a getter is
        # justified in this case.
        if not only_first:
            return self.metadata[label] if label in self.metadata else []
        else:
            return self.metadata[label][0] if self.metadata[label] else None


class FrameshiftNode(object):
    # This class helps us generate all possible frameshift combinations on a transcript.
    # We represent them as a binary tree, where a path from the root to a leaves is a valid
    # frameshift combination. The first level of the tree corresponds to the first frameshift
    # appearing on the transcript, the second level to the second frameshift and so on.
    # The path from the root to a leaf is a sequence of steps each taken either to the left or
    # the right, where a left step means frameshift is "on" and a right step means "off".
    # As homozygous frameshifts can't be "off" their right child is always a None leaf.
    # TODO: explain how the path enumerator function gives you the legal combinations in
    # such an order that you can use it to your advantage when pruning equivalent combinations
    # thanks to early stop codons.

    def __init__(self, fs, next=None, genotype='het'):
        assert genotype in ('hom', 'het'), "Frameshift genotype must be either 'hom' or 'het'"
        self.fs = fs
        self.on = next  # next FrameshiftNode on transcript, None if this is the last one
        self.off = None if genotype=='hom' else next  # if FS is 'hom' it can't be "off"

    def traverse(self):
        tron = self.on.traverse() if self.on is not None else []
        troff = self.off.traverse() if self.off is not None else []

        if self.on is None and self.off is None:
            return [[]]

        myresult = []

        for a in tron:
            myresult.append([self] + a)
        for a in troff:
            myresult.append(a)
        #print 'traversed ', self.fs, myresult
        return myresult

    def __repr__(self):
        return str(self.fs)


class Score(object):
    """
    The class Score holds a structure for scoring a prediction. A score is composed of 3 attributes: the pure score,
    the affinity and the associated rank.
    """
    def __init__(self, method, allele, score, affinity, rank):
        self.method = method
        self.allele = allele
        self.score = score
        self.affinity = affinity
        self.rank = rank


class fred2_attrgetter:
    def __init__(self, attr, *attrs):
        if not attrs:
            if not isinstance(attr, str):
                raise TypeError('attribute name must be a string')
            names = attr.split('.')
            def func(obj):
                for name in names:
                    if name == 'seq':
                        if isinstance(obj, AASequence):
                            obj = getattr(obj.seq, '__str__')()
                        else:
                            obj = getattr(obj, '__str__')()
                        return obj
                    else:
                        obj = getattr(obj, name)
                        return obj
            self._call = func
        else:
            getters = tuple(map(fred2_attrgetter, (attr,) + attrs))
            def func(obj):
                return tuple(getter(obj) for getter in getters)
            self._call = func

    def __call__(self, obj):
        return self._call(obj)

    # for benchmarking
    # from fred2.Core.Base import AASequence
    # from fred2.Core.Base import fred2_attrgetter
    # from fred2 import Core
    # fag = fred2_attrgetter('seq')
    # test = [AASequence('SYFPEITHI','id1','name1','desc1'),AASequence('SYFPEITHI','id2','name2','desc2'),AASequence('SYFPEITHI','id3','name3','desc3'),AASequence('IHTIEPFYS','id1','name1','desc1')]
    # fag(test[0]) ...


class AASequence(SeqRecord):
    def __init__(self, seq, id='<unknown id>', name='<unknown name>', description='<unknown description>',
                 tid=None, pid=None, dbxrefs=None, features=None, annotations=None, letter_annotations=None):
        """
        Forces seq into a Seq object, removing whitespace and other special characters
        :type self: SeqRecord
        """
        aas = Seq(re.sub(r'(\s+|[^\w])', '', seq.upper()), IUPAC.protein) if not isinstance(seq, Seq) else seq
        # TODO warning when seq has IUPAC.protein characters or Bio.Seq is used
        SeqRecord.__init__(self, aas, id, name, description, dbxrefs, features, annotations, letter_annotations)
        self.variants = dict()  # variants
        self.scores = list()  # actually dict(dict())
        self.tid = tid
        self.pid = pid

    def add_snv(self, variant):
        try:
            ori, sub = re.compile(r"\d+").split(variant.coding[self.tid].aa_mutation_syntax.strip("p."))
            pos = int(re.findall(r"\d+", variant.coding[self.tid].aa_mutation_syntax)[0])
            if self.seq[pos-1] == ori:
                if pos-1 in self.variants:
                    raise AttributeError('Already a Variant at this position.')  # TODO handle better & more generic (in/del!)
                self.seq = self.seq[:pos-1] + sub + self.seq[pos:]
                self.description += str(variant)
                self.variants[pos-1] = variant
            else:
                raise AttributeError('No position/aminoacid match, coding is off.')
        except AttributeError as ae:
            raise AttributeError('No coding for given variant.')

    def unfold(self, length):
        if length > len(self.seq):
            raise IndexError('Sequence too short to unfold with that length.')
        ret = list()
        anchors = list()
        for pos, var in self.variants.items():
            if var.specific:
                anchors.append(pos)
        for a in anchors:
            sled_heads = range(max(0, a-(length-1)), min(a+1, (len(self.seq)-length)))
            for sled in sled_heads:
                pep = AASequence(self.seq[sled:sled+length], 'fred2|' + self.pid, "", "peptide around " + str(a))
                for vp in self.variants:
                    if vp < sled+length and vp >= sled:
                        pep.variants[vp-sled] = self.variants[vp]
                ret.append(pep)  # TODO paranoia checks
                # TODO add resp. variants to each fold
        return ret

    # for benchmarking
    #def randomword(length):
    #    return ''.join(random.choice(string.uppercase) for i in range(length))
    # ids = ['gene1','gene2','gene3']
    # tests = [AASequence(randomword(9),i) for i in ids*5]


