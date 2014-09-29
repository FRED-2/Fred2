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
        This class allows automatic registration of new plugins.
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


class ACleavageSitePrediction(object):
    __metaclass__ = APluginRegister

    @abc.abstractproperty
    def name(self):
        """
        Returns the name of the predictor

        :return:
        """
        raise NotImplementedError

    @abc.abstractproperty
    def supportedLength(self):
        """
        Returns the supported lengths of the predictor

        :return: list(int) - Supported peptide length
        """
        raise  NotImplementedError


    @abc.abstractproperty
    def cleavagePos(self):
        """
        parameter specifying the position of aa (within the prediction window) after which the sequence is cleaved
        (starting from 1)

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


class ACleavageFragmentPrediction(object):
    __metaclass__ = APluginRegister

    @abc.abstractproperty
    def name(self):
        """
        Returns the name of the predictor

        :return:
        """
        raise NotImplementedError

    @abc.abstractproperty
    def supportedLength(self):
        """
        Returns the supported lengths of the predictor

        :return: list(int) - Supported peptide length
        """
        raise  NotImplementedError


    @abc.abstractproperty
    def cleavagePos(self):
        """
        parameter specifying the position of aa (within the prediction window) after which the sequence is cleaved

        :return:
        """
        raise NotImplementedError

    @abc.abstractmethod
    def predict(self, _aa_seq,  **kwargs):
        """
        Predicts the probability that the fragment can be produced by the proteasom

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
        Converts alleles into the interal allele representation of the predictor
        and returns a string representation

        :param list(Allele) alleles: The alleles for which the internal predictor
         representation is needed
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


class ASVM(object):
    """
        Base class for SVM prediction tools
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def encode(self, peptides):
        """
        Returns the feature encoding for peptides

        :param List(Peptide)/Peptide peptides: List oder Peptide object
        :return: list(Object) -- Feature encoding of the Peptide objects
        """
        raise NotImplementedError


class AExternal(object):
    """
     Base class for external tools
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def externalPath(self):
        """
        defines the external execution path ?
        :return:
        """

    @abc.abstractmethod
    def parse_external_result(self, _file):
        """
        Parses external NetMHC results and returns a Result object

        :param str _file: The file path or the external prediction results
        :return: Result - Returns a Result object
        """
        raise NotImplementedError

class ATAPPrediction(object):
    __metaclass__ = APluginRegister

    @abc.abstractproperty
    def name(self):
        """
        Returns the name of the predictor

        :return:
        """
        raise NotImplementedError

    @abc.abstractproperty
    def supportedLength(self):
        """
        Returns the supported lengths of the predictor

        :return: list(int) - Supported peptide length
        """
        raise NotImplementedError

    @abc.abstractmethod
    def predict(self, peptides,  **kwargs):
        """
        Predicts the TAP affinity for the given sequences

        :param list(Peptide)/Peptide: Peptides for which TAP affinity should be predicted
        :return: Returns a TAPResult object
        """
        raise NotImplementedError