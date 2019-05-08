# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: TAPPRediction
   :synopsis: Base class for TAP prediction methods
.. moduleauthor:: schubert
"""
from Fred2.Core.Base import ATAPPrediction
from Fred2.TAPPrediction.SVM import *
from Fred2.TAPPrediction.PSSM import *

try:
    from fred_plugin import *
except ImportError:
    pass


class TAPPredictorFactory(object):
    class __metaclass__(type):
        def __init__(cls, name, bases, nmspc):
            type.__init__(cls, name, bases, nmspc)

        def __call__(self, _predictor, *args, **kwargs):
            '''
            If a third person wants to write a new Epitope Predictior. He/She has to name the file fred_plugin and
            inherit from AEpitopePrediction. That's it nothing more.
            '''

            version = str(kwargs["version"]).lower() if "version" in kwargs else None
            try:
                return ATAPPrediction[str(_predictor.lower()), version](*args)
            except KeyError as e:
                if version is None:
                    raise ValueError("Predictor %s is not known. Please verify that such an Predictor is "%_predictor +
                                "supported by FRED2 and inherits ATAPPrediction.")
                else:
                    raise ValueError("Predictor %s version %s is not known. Please verify that such an Predictor is "%(_predictor, version) +
                                "supported by FRED2 and inherits ATAPPrediction.")

    @staticmethod
    def available_methods():
        """
        Returns a dictionary of available TAP predictors and the supported versions

        :return: dict(str, list(str)) - A dictionary of TAP predictors represented as string and supported versions
        """
        return {k:sorted(versions.keys()) for k,versions in ATAPPrediction.registry.items()}