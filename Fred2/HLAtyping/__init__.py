from Fred2.Core.Base import AHLATyping
from Fred2.HLAtyping.External import *

try:
    from fred_plugin import *
except ImportError:
    pass


class HLATypingFactory(object):
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
                return AHLATyping[str(_predictor).lower(), version](*args)
            except KeyError as e:
                if version is None:
                    raise ValueError("Predictor %s is not known. Please verify that such an Predictor is "%_predictor +
                                "supported by FRED2 and inherits AHLATyping.")
                else:
                    raise ValueError("Predictor %s version %s is not known. Please verify that such an Predictor is "%(_predictor, version) +
                                "supported by FRED2 and inherits AHLATyping.")

    @staticmethod
    def available_methods():
        """
        Returns a dictionary of available epitope predictors and the supported versions

        :return: dict(str, list(str)) - A dictionary of epitope predictors represented as string and a list of versions
        """
        return {k:sorted(versions.keys()) for k,versions in AHLATyping.registry.items() }