__author__ = 'schubert'

from Prediction import *
from Prediction import EpitopePredictorFactory


dummy = EpitopePredictorFactory("DummyEpitopePredictor")
#print dummy.__class__
#print dir(dummy)

#print dummy.predict(None)