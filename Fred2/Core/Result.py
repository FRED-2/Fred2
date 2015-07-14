# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.AResult
   :synopsis: Contains relevant classes describing results of predictions.
.. moduleauthor:: schubert,

"""
__author__ = 'schubert'

import abc
import numpy
import pandas


class AResult(pandas.DataFrame):
    """
        A AResult object is a DataFrame with with multi-indexing.

        This class is used as interface and cann be extended with custom short-cuts for the sometimes often tedious
        calls in pandas
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def filter_result(self, expression):
        """
        filter result based on a list of expressions

        :param expression: list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator, threshold)
        :return: AResult
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def merge_results(self, others):
        """
        Merges results of the same type and returns a merged result

        :param AResult other:
        :return:
        """
        raise NotImplementedError()


class EpitopePredictionResult(AResult):
    """
        A AResult object is a DataFrame with with multi-indexing, where column Id are the prediction model (i.e HLA allele
        for epitope prediction), row ID the target of the prediction (i.e. peptide) and the second row ID the predictor
        (i.e BIMAS)

        Epitope prediction result
    """

    def filter_result(self, expressions):
        """
        filters a result data frame based on a specified expression consisting
        of a list of triple with (method_name, comparator, threshold). The expression is applied to each row. If any of the
        columns full fill the criteria the row remains.

        :param list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator, threshold)
        :return: filtered result object
        """
        #builde logical expression
        masks = numpy.logical_and(*map(list,[comp(self.loc[(slice(None), method), :], thr).any(axis=1)
                                             for method, comp, thr in expressions]))

        #apply to all rows
        idx = [f for f in masks
                    for _ in xrange(len(self.index.levels[1]))]
        return EpitopePredictionResult(self.loc[idx, :])

    def merge_results(self, others):
        df = self.copy(deep=False)

        if isinstance(others, EpitopePredictionResult):
            others = [others]

        for i in xrange(len(others)):
            df1a, df2a = df.align(others[i])
            zero1 = df1a == 0
            zero2 = df2a == 0
            df1a = df1a.fillna(0)
            df2a = df2a.fillna(0)
            df = df1a+df2a
            true_zero = zero1 | zero2
            false_zero = df == 0
            zero = true_zero & false_zero
            print zero
            nans = ~true_zero & false_zero
            df[zero] = 0
            df[nans] = numpy.NaN
        return EpitopePredictionResult(df)


class CleavageSitePredictionResult(AResult):
    """
        Epitope prediction result
    """

    def filter_result(self, expressions):
        return None

    def merge_results(self, others):
        if isinstance(others, CleavageSitePredictionResult):
            others = [others]

        df = pandas.concat([self]+others, axis=1)
        return df.T.drop_duplicates().T


class CleavageFragmentPredictionResult(AResult):
    """
        Epitope prediction result
    """

    def filter_result(self, expressions):
        """
        filters a result data frame based on a specified expression consisting
        of a list of triple with (method_name, comparator, threshold). The expression is applied to each row. If any of the
        columns full fill the criteria the row remains.

        :param list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator, threshold)
        :return: filtered result object
        """
        #builde logical expression
        masks = numpy.logical_and(*map(list,[comp(self.loc[(slice(None), method), :], thr).any(axis=1)
                                             for method, comp, thr in expressions]))

        #apply to all rows
        idx = [f for f in masks
                    for _ in xrange(len(self.index.levels[1]))]
        return CleavageFragmentPredictionResult(self.loc[idx, :])

    def merge_results(self, others):
        df = self.copy(deep=False)

        for i in xrange(len(others)):
            df1a, df2a = df.align(others[i])
            zero1 = df1a == 0
            zero2 = df2a == 0
            df1a = df1a.fillna(0)
            df2a = df2a.fillna(0)
            df = df1a+df2a
            true_zero = zero1 | zero2
            false_zero = df == 0
            zero = true_zero & false_zero
            print zero
            nans = ~true_zero & false_zero
            df[zero] = 0
            df[nans] = numpy.NaN
        return CleavageFragmentPredictionResult(df)


class TAPPredictionResult(AResult):
    """
        Epitope prediction result
    """

    def filter_result(self, expressions):
        """
        filters a result data frame based on a specified expression consisting
        of a list of triple with (method_name, comparator, threshold). The expression is applied to each row. If any of the
        columns full fill the criteria the row remains.

        :param list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator, threshold)
        :return: filtered result object
        """
        #builde logical expression
        masks = numpy.logical_and(*map(list,[comp(self.loc[(slice(None), method), :], thr).any(axis=1)
                                             for method, comp, thr in expressions]))

        #apply to all rows
        idx = [f for f in masks
                    for _ in xrange(len(self.index.levels[1]))]
        return TAPPredictionResult(self.loc[idx, :])

    def merge_results(self, others):
        df = self.copy(deep=False)

        for i in xrange(len(others)):
            df1a, df2a = df.align(others[i])
            zero1 = df1a == 0
            zero2 = df2a == 0
            df1a = df1a.fillna(0)
            df2a = df2a.fillna(0)
            df = df1a+df2a
            true_zero = zero1 | zero2
            false_zero = df == 0
            zero = true_zero & false_zero
            print zero
            nans = ~true_zero & false_zero
            df[zero] = 0
            df[nans] = numpy.NaN
        return TAPPredictionResult(df)
