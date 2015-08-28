# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: CleavagePrediction
   :synopsis: Factory classes for cleavage site and fragment prediction.
              This is the entry point to all cleavage prediction methods.
.. moduleauthor:: schubert

"""

from Fred2.Core.Base import ACleavageSitePrediction, ACleavageFragmentPrediction
from Fred2.CleavagePrediction.PSSM import *
from Fred2.CleavagePrediction.External import *

try:
    from fred_plugin import *
except ImportError:
    pass


class CleavageSitePredictorFactory(object):
    class __metaclass__(type):
        def __init__(cls, name, bases, nmspc):
            type.__init__(cls, name, bases, nmspc)

        def __call__(self, _predictor, *args, **kwargs):
            '''
            just as I think it works.....

            If a third person wants to write a new Cleavage Site Predictor. One has to name the file fred_plugin and
            inherit from ACleavagePrediction. That's it nothing more.
            '''

            try:
                return ACleavageSitePrediction.registry[_predictor](*args)
            except KeyError:
                raise ValueError("Predictor %s is not known. Please verify that such an predictor is "%_predictor +
                                 "supported by FRED2 and inherits ACleavageSitePrediction.")

    @staticmethod
    def available_methods():
        """
        Returns a list of available cleavage site predictors

        :return: list of cleavage site predictor represented as string
        """
        return sorted([k for k in ACleavageSitePrediction.registry.iterkeys()
                       if k not in ["APSSMCleavageSitePredictor", "ACleavageSitePrediction"]])


class CleavageFragmentPredictorFactory(object):
    class __metaclass__(type):
        def __init__(cls, name, bases, nmspc):
            type.__init__(cls, name, bases, nmspc)

        def __call__(self, _predictor, *args, **kwargs):
            '''
            just as I think it works.....

            If a third person wants to write a new Cleavage Fragment Predictor. One has to name the file fred_plugin and
            inherit from ACleavagePrediction. That's it nothing more.
            '''

            try:
                return ACleavageFragmentPrediction.registry[_predictor](*args)
            except KeyError:
                raise ValueError("Predictor %s is not known. Please verify that such an predictor is "%_predictor +
                                 "supported by FRED2 and inherits ACleavageFragmentPrediction.")

    @staticmethod
    def available_methods():
        """
        Returns a list of available cleavage fragment predictors

        :return: list of cleavage fragment represented as string
        """
        return sorted([k for k in ACleavageFragmentPrediction.registry.iterkeys()
                if k not in ["APSSMCleavageFragmentPredictor", "ACleavageFragmentPrediction"]])