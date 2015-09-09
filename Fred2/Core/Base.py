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
        if not inspect.isabstract(cls):
            cls.registry.setdefault(str(cls().name).lower(), {}).update({str(cls().version).lower():cls})

    def __getitem__(cls, args):
        name, version = args
        if version is None:
            return cls.registry[name][max(cls.registry[name].keys())]
        return cls.registry[name][version]

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
    def version(self):
        """
        parameter specifying the version of the prediction method

        """
        raise NotImplementedError

    @abc.abstractproperty
    def supportedLength(self):
        """
        Returns the supported lengths of the predictor

        :return: list(int) - Supported peptide length
        """
        raise NotImplementedError

    @abc.abstractproperty
    def cleavagePos(self):
        """
        parameter specifying the position of aa (within the prediction window) after which the sequence is cleaved
        (starting from 1)

        :return:
        """
        raise NotImplementedError


    @abc.abstractmethod
    def predict(self, _aa_seq, **kwargs):
        """
        Predicts the proteom cleavage site of the given sequences

        :param Bio.Seq _aa_seq: The sequence to be cleaved (must be an instance of Bio.Seq
        :return: Returns a AResult object for the specified Bio.Seq
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
    def version(self):
        """
        parameter specifying the version of the prediction method

        """
        raise NotImplementedError

    @abc.abstractproperty
    def supportedLength(self):
        """
        Returns the supported lengths of the predictor

        :return: list(int) - Supported peptide length
        """
        raise NotImplementedError


    @abc.abstractproperty
    def cleavagePos(self):
        """
        parameter specifying the position of aa (within the prediction window) after which the sequence is cleaved

        :return:
        """
        raise NotImplementedError

    @abc.abstractmethod
    def predict(self, _aa_seq, **kwargs):
        """
        Predicts the probability that the fragment can be produced by the proteasom

        :param Bio.Seq _aa_seq: The sequence to be cleaved (must be an instance of Bio.Seq
        :return: Returns a AResult object for the specified Bio.Seq
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
    def version(cls):
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
        :return: Returns a AResult object for the specified Peptides and Alleles
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
    def command(self):
        """
        defines the external execution path ?
        :return:
        """

    @abc.abstractmethod
    def parse_external_result(self, _file):
        """
        Parses external NetMHC results and returns a AResult object

        :param str _file: The file path or the external prediction results
        :return: AResult - Returns a AResult object
        """
        raise NotImplementedError

    def is_in_path(self):
        """
        checks whether the specified execution command can be found in PATH

        :return: bool - Whether or not command could be found in PATH
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

        :param (str) path: - optional specification of executable path if deviant from self.__command
        :return: str - The external version of the tool
        """
        exe = self.command.split()[0] if path is None else path
        try:
            p = subprocess.Popen(exe + ' --version', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait() #block the rest
            stdo, stde = p.communicate()
            stdr = p.returncode
            if stdr > 0:
                raise RuntimeError("Could not check version of " + exe + " - Please check your installation and FRED2 "
                                                                         "wrapper implementation.")
        except Exception as e:
                raise RuntimeError(e)
        return str(stdo).strip()

    @abc.abstractmethod
    def prepare_peptide_input(self, _peptides, _file):
        """
        Prepares sequence input for external tools
        and writes them to _file in the specific format

        NO return value!

        :param: (list(str)) _peptides: the peptide sequences to write into _file
        :param (File) _file: File handler to input file for external tool
        """
        return NotImplementedError


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
    def version(self):
        """
        parameter specifying the version of the prediction method

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
    def predict(self, peptides, **kwargs):
        """
        Predicts the TAP affinity for the given sequences

        :param list(Peptide)/Peptide: Peptides for which TAP affinity should be predicted
        :return: Returns a TAPResult object
        """
        raise NotImplementedError