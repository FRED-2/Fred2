# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.Base
   :synopsis: This module contains base classes for all other modules.
.. moduleauthor:: schubert, szolek, walzer


https://docs.python.org/3/library/abc.html

"""


import abc
import inspect
import os
import subprocess
from collections import defaultdict
from string import maketrans
import warnings


COMPLEMENT = maketrans('atgcATGC', 'tacgTACG')


class MetadataLogger(object):
    """
    This class provides a simple interface for assigning additional metadata to
    any object in our data model. Examples: storing ANNOVAR columns like depth,
    base count, dbSNP id, quality information for variants, additional prediction information
    for peptides etc. This functionality is not used from core methods of FRED2.

    The saved values are accessed via :meth:`~Fred2.Core.MetadataLogger.log_metadata` and
    :meth:`~Fred2.Core.MetadataLogger.get_metadata`
    
    """
    def __init__(self):
        """
        """
        self.__metadata = defaultdict(list)

    def log_metadata(self, label, value):
        """
        Inserts a new metadata

        :param str label: key for the metadata that will be added
        :param list(object) value: any kindy of additional value that should be kept
        """
        self.__metadata[label].append(value)

    def get_metadata(self, label, only_first=False):
        """
        Getter for the saved metadata with the key :attr:`label`

        :param str label: key for the metadata that is inferred
        :param bool only_first: true if only the the first element of the matadata list is to be returned
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
        if not inspect.isabstract(cls):
            cls.registry.setdefault(str(cls().name).lower(), {}).update({str(cls().version).lower():cls})

    def __getitem__(cls, args):
        name, version = args
        if version is None:
            return cls.registry[name][max(cls.registry[name].keys())]
        return cls.registry[name][version]

    def __iter__(cls):
        return iter(list(cls.registry.values()))

    def __str__(cls):
        if cls in cls.registry:
            return cls.__name__
        return cls.__name__ + ": " + ", ".join([sc.__name__ for sc in cls])


class ACleavageSitePrediction(object, metaclass=APluginRegister):
    @abc.abstractproperty
    def name(self):
        """
        The name of the predictor
        """
        raise NotImplementedError

    @abc.abstractproperty
    def version(self):
        """
        Parameter specifying the version of the prediction method
        """
        raise NotImplementedError

    @abc.abstractproperty
    def supportedLength(self):
        """
        The supported lengths of the predictor
        """
        raise NotImplementedError

    @abc.abstractproperty
    def cleavagePos(self):
        """
        Parameter specifying the position of aa (within the prediction window) after which the sequence is cleaved
        (starting from 1)
        """
        raise NotImplementedError


    @abc.abstractmethod
    def predict(self, aa_seq, **kwargs):
        """
        Predicts the proteasomal cleavage site of the given sequences

        :param aa_seq: The sequence to be cleaved
        :type aa_seq: :class:`~Fred2.Core.Peptide.Peptide` or :class:`~Fred2.Core.Protein.Protein`
        :return: Returns a :class:`~Fred2.Core.Result.AResult` object for the specified Bio.Seq
        :rtype: :class:`~Fred2.Core.Result.AResult`
        """
        raise NotImplementedError


class ACleavageFragmentPrediction(object, metaclass=APluginRegister):
    @abc.abstractproperty
    def name(self):
        """
        The name of the predictor
        """
        raise NotImplementedError

    @abc.abstractproperty
    def version(self):
        """
        Parameter specifying the version of the prediction method
        """
        raise NotImplementedError

    @abc.abstractproperty
    def supportedLength(self):
        """
        The supported lengths of the predictor
        """
        raise NotImplementedError


    @abc.abstractproperty
    def cleavagePos(self):
        """
        Parameter specifying the position of aa (within the prediction window) after which the sequence is cleaved
        """
        raise NotImplementedError

    @abc.abstractmethod
    def predict(self, aa_seq, **kwargs):
        """
        Predicts the probability that the fragment can be produced by the proteasom

        :param aa_seq: The sequence to be cleaved
        :type aa_seq: :class:`~Fred2.Core.Peptide.Peptide`
        :return: Returns a :class:`~Fred2.Core.Result.AResult` object for the specified Bio.Seq
        :rtype: :class:`~Fred2.Core.Result.AResult`
        """
        raise NotImplementedError


class AEpitopePrediction(object, metaclass=APluginRegister):
    @abc.abstractproperty
    def name(self):
        """
        The name of the predictor
        """
        raise NotImplementedError

    @abc.abstractproperty
    def version(cls):
        """The version of the predictor"""
        raise NotImplementedError

    @abc.abstractproperty
    def supportedAlleles(self):
        """
        A list of valid allele models
        """
        raise NotImplementedError

    @abc.abstractproperty
    def supportedLength(self):
        """
        A list of supported peptide lengths
        """
        raise NotImplementedError

    @abc.abstractmethod
    def convert_alleles(self, alleles):
        """
        Converts alleles into the internal allele representation of the predictor
        and returns a string representation

        :param alleles: The alleles for which the internal predictor representation is needed
        :type alleles: list(:class:`~Fred2.Core.Allele.Allele`)
        :return: Returns a string representation of the input alleles
        :rtype: list(str)
        """
        raise NotImplementedError

    @abc.abstractmethod
    def predict(self, peptides, alleles=None, **kwargs):
        """
        Predicts the binding affinity for a given peptide or peptide lists for a given list of alleles.
        If alleles is not given, predictions for all valid alleles of the predictor is performed. If, however,
        a list of alleles is given, predictions for the valid allele subset is performed.

        :param peptides: The peptide objects for which predictions should be performed
        :type peptides: :class:`~Fred2.Core.Peptide.Peptide` or list(:class:`~Fred2.Core.Peptide.Peptide`)
        :param alleles: An :class:`~Fred2.Core.Allele.Allele` or list of :class:`~Fred2.Core.Allele.Allele` for which
                        prediction models should be used
        :type alleles: :class:`~Fred2.Core.Allele.Allele`/list(:class:`~Fred2.Core.Allele.Allele`)
        :return: Returns a :class:`~Fred2.Core.Result.AResult` object for the specified
                 :class:`~Fred2.Core.Peptide.Peptide` and :class:`~Fred2.Core.Allele.Allele`
        :rtype: :class:`~Fred2.Core.Result.AResult`
        """
        raise NotImplementedError


class ASVM(object, metaclass=abc.ABCMeta):
    """
        Base class for SVM prediction tools
    """

    @abc.abstractmethod
    def encode(self, peptides):
        """
        Returns the feature encoding for peptides

        :param peptides: List of or a single :class:`~Fred2.Core.Peptide.Peptide` object
        :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`)/:class:`~Fred2.Core.Peptide.Peptide`
        :return: Feature encoding of the Peptide objects
        :rtype: list(Object)
        """
        raise NotImplementedError


class AExternal(object, metaclass=abc.ABCMeta):
    """
     Base class for external tools
    """

    @abc.abstractproperty
    def command(self):
        """
        Defines the commandline call for external tool
        """
        raise NotImplementedError

    @abc.abstractmethod
    def parse_external_result(self, file):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :return: A dictionary containing the prediction results
        :rtype: dict
        """
        raise NotImplementedError

    def is_in_path(self):
        """
        Checks whether the specified execution command can be found in PATH

        :return: Whether or not command could be found in PATH
        :rtype: bool
        """
        exe = self.command.split()[0]
        for try_path in os.environ["PATH"].split(os.pathsep):
            try_path = try_path.strip('"')
            exe_try = os.path.join(try_path, exe).strip()
            if os.path.isfile(exe_try) and os.access(exe_try, os.X_OK):
                return True
        return False

    @abc.abstractmethod
    def get_external_version(self, path=None):
        """
        Returns the external version of the tool by executing
        >{command} --version

        might be dependent on the method and has to be overwritten
        therefore it is declared abstract to enforce the user to
        overwrite the method. The function in the base class can be called
        with super()

        :param str path: - Optional specification of executable path if deviant from self.__command
        :return: The external version of the tool or None if tool does not support versioning
        :rtype: str
        """
        exe = self.command.split()[0] if path is None else path
        try:
            p = subprocess.Popen(exe + ' --version', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #p.wait() #block the rest
            stdo, stde = p.communicate()
            stdr = p.returncode
            if stdr > 0:
                raise RuntimeError("Could not check version of " + exe + " - Please check your installation and FRED2 "
                                                                         "wrapper implementation.")
        except Exception as e:
                raise RuntimeError(e)
        return str(stdo).strip()


class ATAPPrediction(object, metaclass=APluginRegister):
    @abc.abstractproperty
    def name(self):
        """
        The name of the predictor
        """
        raise NotImplementedError

    @abc.abstractproperty
    def version(self):
        """
        Parameter specifying the version of the prediction method
        """
        raise NotImplementedError

    @abc.abstractproperty
    def supportedLength(self):
        """
        The supported lengths of the predictor
        """
        raise NotImplementedError


    @abc.abstractmethod
    def predict(self, peptides, **kwargs):
        """
        Predicts the TAP affinity for the given sequences

        :param peptides: :class:`~Fred2.Core.Peptide.Peptide` for which TAP affinity should be predicted
        :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`)/:class:`~Fred2.Core.Peptide.Peptide`
        :return: Returns a :class:`~Fred2.Core.Result.TAPResult` object
        :rtype: :class:`~Fred2.Core.Result.TAPResult`
        """
        raise NotImplementedError


class AHLATyping(object, metaclass=APluginRegister):
    @abc.abstractproperty
    def name(self):
        """
        The name of the predictor
        """
        raise NotImplementedError

    @abc.abstractproperty
    def version(self):
        """
        Parameter specifying the version of the prediction method
        """
        raise NotImplementedError

    @abc.abstractmethod
    def predict(self, ngsFile, output, **kwargs):
        """
        Prediction method for inferring the HLA typing

        :param str ngsFile: The path to the input file containing the NGS reads
        :param str output: The path to the output file or directory
        :return: A list of HLA alleles representing the genotype predicted by the algorithm
        :rtype: list(:class:`~Fred2.Core.Allele.Allele`)
        """
        raise NotImplementedError


def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""
    def new_func(*args, **kwargs):
        warnings.simplefilter('default')  #this will render these deprecation warnings visible to everyone (default is switched off in python >=2.7)
        warnings.warn("Call to deprecated function {n} of {f}.".format(n=func.__name__, f=func.__doc__),
                      category=DeprecationWarning)
        return func(*args, **kwargs)
    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func