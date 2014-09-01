# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'schubert'

import abc

from Fred2.Prediction.NetMHC import *
from Fred2.Prediction.PSSM import *
from Fred2.Core.Base import PluginRegister



class AEpitopePrediction(object):
    __metaclass__ = PluginRegister

    @abc.abstractproperty
    def alleleModels(self):
        """
        Returns a list of valid allele models

        :return: List of allele names for which the predictor provides models
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

        
class EpitopePredictorFactory(object):
    class __metaclass__(type):
        def __init__(cls, name, bases, nmspc):
            type.__init__(cls, name, bases, nmspc)

        def __call__(self, _predictor, *args, **kwargs):
            '''
            just as I think it works.....

            If a third person wants to write a new Epitope Predictior. He/She has to name the file fred_plugin and
            inherit from AEpitopePrediction. That's it nothing more.
            '''
            try:
                from fred_plugin import *
            except ImportError:
                pass

            try:
                return AEpitopePrediction.registry[_predictor](*args)
            except KeyError:
                raise ValueError("Predictor %s is not known. Please verify that such an Predictor is "%_predictor +
                                 "supported by FRED2 and inherits AEpitopePredictor.")
