__author__ = 'schubert'

from EpitopePrediction import *
from EpitopePrediction import EpitopePredictorFactory


dummy = EpitopePredictorFactory("DummyEpitopePredictor")
#print dummy.__class__
#print dir(dummy)

#print dummy.predict(None)