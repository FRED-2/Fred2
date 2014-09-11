# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'schubert'

from Fred2.Core.Base import AEpitopePrediction
from Fred2.EpitopePrediction.ANN import *
from Fred2.EpitopePrediction.PSSM import *
from Fred2.EpitopePrediction.SVM import *

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

    @staticmethod
    def available_methods(self):
        """
        Returns a list of available epitope predictors

        :return: list of epitope predictors represented as string
        """
        return AEpitopePrediction.registry.keys()