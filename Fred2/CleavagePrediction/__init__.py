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
            """
            If a third person wants to write a new Cleavage Site Predictor. One has to name the file fred_plugin and
            inherit from ACleavagePrediction. That's it nothing more.
            """

            version = str(kwargs["version"]).lower() if "version" in kwargs else None
            try:
                return ACleavageSitePrediction[str(_predictor).lower(), version](*args)
            except KeyError as e:
                if version is None:
                    raise ValueError("Predictor %s is not known. Please verify that such an Predictor is "%_predictor +
                                "supported by FRED2 and inherits ACleavageSitePrediction.")
                else:
                    raise ValueError("Predictor %s version %s is not known. Please verify that such an Predictor is "%(_predictor, version) +
                                "supported by FRED2 and inherits ACleavageSitePrediction.")

    @staticmethod
    def available_methods():
        """
        Returns a list of available cleavage site predictors

        :return: dict(str,list(int)) - dict of cleavage site predictor represented as string and the supported versions
        """
        return {k: sorted(versions.keys()) for k, versions in ACleavageSitePrediction.registry.items()}


class CleavageFragmentPredictorFactory(object):
    class __metaclass__(type):
        def __init__(cls, name, bases, nmspc):
            type.__init__(cls, name, bases, nmspc)

        def __call__(self, _predictor, *args, **kwargs):
            """
            If a third person wants to write a new Cleavage Fragment Predictor. One has to name the file fred_plugin and
            inherit from CleavageFragmentPredictorFactory. That's it nothing more.
            """

            version = str(kwargs["version"]).lower() if "version" in kwargs else None
            try:
                return ACleavageFragmentPrediction[str(_predictor).lower(), version](*args)
            except KeyError as e:
                if version is None:
                    raise ValueError("Predictor %s is not known. Please verify that such an Predictor is "%_predictor +
                                "supported by FRED2 and inherits ACleavageFragmentPrediction.")
                else:
                    raise ValueError("Predictor %s version %s is not known. Please verify that such an Predictor is "%(_predictor, version) +
                                "supported by FRED2 and inherits ACleavageFragmentPrediction.")

    @staticmethod
    def available_methods():
        """
        Returns a list of available cleavage site predictors

        :return: dict(str,list(str)) - dict of cleavage site predictor represented as string and the supported versions
        """
        return {k: sorted(versions.keys()) for k, versions in ACleavageFragmentPrediction.registry.items()}