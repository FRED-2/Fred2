__author__ = 'schubert'

'''
This is just a dummy
'''
from Prediction import AEpitopePrediction


class DummyEpitopePredictor(AEpitopePrediction):

    __alleles = ["HLA-A*01:01"]
    __name__ = "DummyPredictor"

    @property
    def supportedAlleles(self):
        return self.__alleles

    def predict(self, peptides, alleles=None, **kwargs):
        return True

    def convert_alleles(self, alleles):
        """
        Converts alleles into the interal allele representation of the predictor and returns a string representation

        :param list(Allele) alleles: The alleles for which the internal predictor representation is needed
        :return: Returns a string representation of the input alleles
        """
        return self.__alleles