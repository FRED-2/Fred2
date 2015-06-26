# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.Result
   :synopsis: Contains relevant classes describing results of predictions.
.. moduleauthor:: schubert,

"""
__author__ = 'schubert'

import abc
from pandas import DataFrame


class Result(DataFrame):
    """
        A Result object is a DataFrame with with multi-indexing.

        This class is used as interface and cann be extended with custom short-cuts for the sometimes often tedious
        calls in pandas
    """


class EpitopePredictionResult(Result):
    """
        A Result object is a DataFrame with with multi-indexing, where column Id are the prediction model (i.e HLA allele
        for epitope prediction), row ID the target of the prediction (i.e. peptide) and the second row ID the predictor
        (i.e BIMAS)

        Epitope prediction result
    """

    @staticmethod
    def merge_result(self, _epi_pred_result):
        """
        Takes a epitope result object and merges the two and returns the merge
        result object

        :param: (EpitopePredictionResult) _epi_pred_result: the result object to be merged
        :return: EpitopePredictionResult - the merged EpitopePredictionResult object
        :raise: ValueError - if _epi_pred_result is not of type EpitopePredictionResult
        """
        if not isinstance(_epi_pred_result, EpitopePredictionResult):
            raise ValueError("Input is not of type EpitopePredictionResult.")

        df1a, df2a = self.align(_epi_pred_result)
        return EpitopePredictionResult(df1a+df2a)

    @staticmethod
    def filter_result(self, _filter_criteria):
        """
        Takes a list triples (Prediction_Method,Threshold,Sign) where Prediction_Method is the prediction methods
        name to filter for, Threshold the value an epitope as to exceed/
        """

class CleavageSitePredictionResult(Result):
    """
        Epitope prediction result
    """

class CleavageFragmentPredictionResult(Result):
    """
        Epitope prediction result
    """

class TAPPredictionResult(Result):
    """
        Epitope prediction result
    """