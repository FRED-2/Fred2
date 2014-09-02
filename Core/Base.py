# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Base
   :synopsis: This module contains base classes for all other modules.
.. moduleauthor:: schubert, szolek, walzer

"""
__author__ = 'szolek', 'walzer'


import abc
from collections import defaultdict
from string import maketrans



COMPLEMENT = maketrans('atgcATGC', 'tacgTACG')


class MetadataLogger(object):
    """
    This class provides a simple interface for assigning additional metadata to
    any object in our data model. Examples: storing ANNOVAR columns like depth,
    base count, dbSNP id, quality information for variants, prediction scores 
    for peptides etc. Normally used by custom toolbox functions and importers.

    The saved values are accessed via :meth:`~Fred2.Core.Base.log_metadata` and
    :meth:`~Fred2.Core.Base.get_metadata`
    
    """
    def __init__(self):
        """
        """
        self.__metadata = defaultdict(list)

    def log_metadata(self, label, value):
        """
        Inserts a new metadata

        :param str label: key for the metadata that will be added
        :param list(object) value: any kindy of additional value that should be
                                   kept
        """
        self.__metadata[label].append(value)

    def get_metadata(self, label, only_first=False):
        """
        Getter for the saved metadata with the key :attr:`label`

        :param str label: key for the metadata that is inferred
        :param bool only_first: true if only the the first element of the 
                                matadata list is to be returned
        """
        # although defaultdict *would* return [] if it didn't find label in 
        # self.metadata, it would come with the side effect of adding label as 
        #a key to the defaultdict, so a getter is justified in this case.
        if not only_first:
            return self.__metadata[label] if label in self.__metadata else []
        else:
            return self.__metadata[label][0] if self.__metadata[label] else None


#Metaclass for Plugins
class APluginRegister(abc.ABCMeta):
    """
        This class is an allows automatic registration of new plugins.
    """

    def __init__(cls, name, bases, nmspc):
        super(APluginRegister, cls).__init__(name, bases, nmspc)

        if not hasattr(cls, 'registry'):
            cls.registry = dict()
        cls.registry[name] = cls

        #dont know if needed
        #del cls.registry[bases.__name__] # Remove base classes

    def __getitem__(cls, item):
        return cls.registry[item]

    def __iter__(cls):
        return iter(cls.registry.values())

    def __str__(cls):
        if cls in cls.registry:
            return cls.__name__
        return cls.__name__ + ": " + ", ".join([sc.__name__ for sc in cls])


class ACleavagePrediction(object):
    __metaclass__ = APluginRegister

    @abc.abstractproperty
    def name(self):
        """
        Returns the name of the predictor

        :return:
        """
        raise NotImplementedError

    @abc.abstractmethod
    def predict(self, _aa_seq,  **kwargs):
        """
        Predicts the proteom cleavage site of the given sequences

        :param Bio.Seq _aa_seq: The sequence to be cleaved (must be an instance of Bio.Seq
        :return: Returns a Result object for the specified Bio.Seq
        """
        raise NotImplementedError

class AEpitopePrediction(object):
    __metaclass__ = APluginRegister

    @abc.abstractproperty
    def name(self):
        """
        Returns the name of the predictor

        :return:
        """
        raise NotImplementedError

    @abc.abstractproperty
    def supportedAlleles(self):
        """
        Returns a list of valid allele models

        :return: List of allele names for which the predictor provides models
        """
        raise NotImplementedError

    @abc.abstractproperty
    def supportedLength(self):
        """
        Returns a list of supported peptide lenghts

        :return: List of supported peptide lengths
        """
        raise NotImplementedError

    @abc.abstractmethod
    def convert_alleles(self, alleles):
        """
        Converts alleles into the interal allele representation of the predictor and returns a string representation

        :param list(Allele) alleles: The alleles for which the internal predictor representation is needed
        :return: Returns a string representation of the input alleles
        """
        raise NotImplementedError

    @abc.abstractmethod
    def predict(self, peptides, alleles=None, **kwargs):
        """
        Predicts the binding affinity for a given peptide or peptide lists for a given list of alleles.
        If alleles is not given, predictions for all valid alleles of the predictor is performed. If, however,
        a list of alleles is given, predictions for the valid allele subset is performed.

        :param Peptide/list(Peptide) peptides: The peptide objects for which predictions should be performed
        :param Allele/list(Allele) alleles: An Allele or list of Alleles for which prediction models should be used
        :return: Returns a Result object for the specified Peptides and Alleles
        """
        raise NotImplementedError

#DEPRECATED
#
# class FrameshiftNode(object):
#     # This class helps us generate all possible frameshift combinations on a transcript.
#     # We represent them as a binary tree, where a path from the root to a leaves is a valid
#     # frameshift combination. The first level of the tree corresponds to the first frameshift
#     # appearing on the transcript, the second level to the second frameshift and so on.
#     # The path from the root to a leaf is a sequence of steps each taken either to the left or
#     # the right, where a left step means frameshift is "on" and a right step means "off".
#     # As homozygous frameshifts can't be "off" their right child is always a None leaf.
#     # TODO: explain how the path enumerator function gives you the legal combinations in
#     # such an order that you can use it to your advantage when pruning equivalent combinations
#     # thanks to early stop codons.
#
#     def __init__(self, fs, next=None, genotype='het'):
#         assert genotype in ('hom', 'het'), "Frameshift genotype must be either 'hom' or 'het'"
#         self.fs = fs
#         self.on = next  # next FrameshiftNode on transcript, None if this is the last one
#         self.off = None if genotype=='hom' else next  # if FS is 'hom' it can't be "off"
#
#     def traverse(self):
#         tron = self.on.traverse() if self.on is not None else []
#         troff = self.off.traverse() if self.off is not None else []
#
#         if self.on is None and self.off is None:
#             return [[]]
#
#         myresult = []
#
#         for a in tron:
#             myresult.append([self] + a)
#         for a in troff:
#             myresult.append(a)
#         #print 'traversed ', self.fs, myresult
#         return myresult
#
#     def __repr__(self):
#         return str(self.fs)
#
#
# class Score(object):
#     """
#     The class Score holds a structure for scoring a prediction. A score is composed of 3 attributes: the pure score,
#     the affinity and the associated rank.
#     """
#     def __init__(self, method, allele, score, affinity, rank):
#         self.method = method
#         self.allele = allele
#         self.score = score
#         self.affinity = affinity
#         self.rank = rank
#
#     def __str__(self):
#         return ','.join([self.method, str(self.allele), str(self.score), str(self.affinity), (str(self.rank) if self.rank else '-')])
#
#
# class fred2_attrgetter:
#     def __init__(self, attr, *attrs):
#         if not attrs:
#             if not isinstance(attr, str):
#                 raise TypeError('attribute name must be a string')
#             names = attr.split('.')
#             def func(obj):
#                 for name in names:
#                     if name == 'seq':
#                         if isinstance(obj, AASequence):
#                             obj = getattr(obj.seq, '__str__')()
#                         else:
#                             obj = getattr(obj, '__str__')()
#                         return obj
#                     else:
#                         obj = getattr(obj, name)
#                         return obj
#             self._call = func
#         else:
#             getters = tuple(map(fred2_attrgetter, (attr,) + attrs))
#             def func(obj):
#                 return tuple(getter(obj) for getter in getters)
#             self._call = func
#
#     def __call__(self, obj):
#         return self._call(obj)
#
#     # for benchmarking
#     # from fred2.Core.Base import AASequence
#     # from fred2.Core.Base import fred2_attrgetter
#     # from fred2 import Core
#     # fag = fred2_attrgetter('seq')
#     # test = [AASequence('SYFPEITHI','id1','name1','desc1'),AASequence('SYFPEITHI','id2','name2','desc2'),AASequence('SYFPEITHI','id3','name3','desc3'),AASequence('IHTIEPFYS','id1','name1','desc1')]
#     # fag(test[0]) ...
#
#
# class AASequence(SeqRecord):
#     def __init__(self, seq, id='<unknown id>', name='<unknown name>', description='<unknown description>',
#                  tid=None, pid=None, dbxrefs=None, features=None, annotations=None, letter_annotations=None):
#         """
#         Forces seq into a Seq object, removing whitespace and other special characters
#         :type self: SeqRecord
#         """
#         aas = Seq(re.sub(r'(\s+|[^\w])', '', seq.upper()), IUPAC.protein) if not isinstance(seq, Seq) else seq
#         # TODO warning when seq has IUPAC.protein characters or Bio.Seq is used
#         SeqRecord.__init__(self, aas, id, name, description, dbxrefs, features, annotations, letter_annotations)
#         self.variants = dict()  # variants
#         self.scores = list()  # actually dict(dict())
#         self.tid = tid
#         self.pid = pid
#
#     # does complicate things, should be in transcript (i.e. a generator using variants (specialized) )
#     def add_snv(self, variant):
#         try:
#             ori, sub = re.compile(r"\d+").split(variant.coding[self.tid].aa_mutation_syntax.strip("p."))
#             pos = int(re.findall(r"\d+", variant.coding[self.tid].aa_mutation_syntax)[0])
#             if self.seq[pos-1] == ori:
#                 if pos-1 in self.variants:
#                     raise AttributeError('Already a Variant at this position.')  # TODO handle better & more generic (in/del!)
#                 self.seq = self.seq[:pos-1] + sub + self.seq[pos:]
#                 self.description += str(variant) + '=' + str(variant.coding[self.tid].aa_mutation_syntax)
#                 self.variants[pos-1] = variant
#             else:
#                 raise AttributeError('No position/aminoacid match, coding is off.')
#         except AttributeError as ae:
#             raise AttributeError('Bad coding for given variant: ' + variant.coding[self.tid].aa_mutation_syntax + ' @ ' + self.tid + '/' + self.pid + '\n' + str(ae))
#
#     def unfold(self, length):
#         if length > len(self.seq):
#             raise IndexError('Sequence too short to unfold with that length.')
#         ret = list()
#         anchors = list()
#         for pos, var in self.variants.items():
#             if var.specific:
#                 anchors.append(pos)
#         for a in anchors:
#             sled_heads = range(max(0, a-(length-1)), min(a+1, (len(self.seq)-length)))
#             for sled in sled_heads:
#                 pep = AASequence(self.seq[sled:sled+length], 'fred2|' + self.pid, "", "peptide around " + str(a))
#                 for vp in self.variants:
#                     if vp < sled+length and vp >= sled:
#                         pep.variants[vp-sled] = self.variants[vp]
#                 ret.append(pep)  # TODO paranoia checks
#                 # TODO add resp. variants to each fold
#         return ret
#
#     # for benchmarking
#     #def randomword(length):
#     #    return ''.join(random.choice(string.uppercase) for i in range(length))
#     # ids = ['gene1','gene2','gene3']
#     # tests = [AASequence(randomword(9),i) for i in ids*5]
#
#
