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

        This class is used as interface and can be extended with custom short-cuts for the sometimes often tedious
        calls in pandas
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def filter_result(self, expression):
        """
        Filter result based on a list of expressions

        :param list((str,comparator,float)) expressions: A list of triples consisting of
                                                        (method_name, comparator, threshold)
        :return: AResult - A new filtered AResult object
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def merge_results(self, others):
        """
        Merges results of the same type and returns a merged result

        :param AResult other: Another AResult object of the same class
        :return: AResult - A new merged AResult object
        """
        raise NotImplementedError()


class EpitopePredictionResult(AResult):
    """
        A AResult object is a DataFrame with with multi-indexing, where column Id are the prediction model (i.e HLA allele
        for epitope prediction), row ID the target of the prediction (i.e. peptide) and the second row ID the predictor
        (e.g. BIMAS)

        Epitope prediction result

        +--------------+-------------+-------------+-------------+-------------+
        | Peptide Obj  | Method Name | Allele1 Obj | Allele2 Obj | Allele3 Obj |
        +==============+=============+=============+=============+=============+
        | Peptide1     | Method 1    |    0.324    |     0.56    |    0.013    |
        +              +-------------+-------------+-------------+-------------+
        |              | Method 2    |     20      |      15     |     23      |
        +--------------+-------------+-------------+-------------+-------------+
        | Peptide2     | Method 1    |    0.50     |     0.36    |    0.98     |
        +              +-------------+-------------+-------------+-------------+
        |              | Method 2    |     26      |      10     |     50      |
        +--------------+-------------+-------------+-------------+-------------+

    """

    def filter_result(self, expressions):
        """
        Filters a result data frame based on a specified expression consisting
        of a list of triple with (method_name, comparator, threshold). The expression is applied to each row. If any of the
        columns full fill the criteria the row remains.

        :param list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator, threshold)
        :return: EpitopePredictionResult - Filtered result object
        """
        if isinstance(expressions, tuple):
            expressions = [expressions]

        #builde logical expression
        masks = map(list, [comp(self.loc[(slice(None), method), :], thr).any(axis=1)
                    for method, comp, thr in expressions])
        if len(masks) > 1:
            masks = numpy.logical_and(*masks)
        else:
            masks = masks[0]

        #apply to all rows
        idx = [f for f in masks
               for _ in xrange(len(self.index.levels[1]))]
        return EpitopePredictionResult(self.loc[idx, :])

    def merge_results(self, others):
        """
        Merges results of type EpitopePredictionResult and returns the merged result

        :param EpitopePredictionResult other: Another EpitopePredictionResult object of the same class
        :return: EpitopePredictionResult - A new merged EpitopePredictionResult object
        """
        df = self.copy(deep=False)

        if type(others) == type(self):
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
            nans = ~true_zero & false_zero
            df[zero] = 0
            df[nans] = numpy.NaN
        return EpitopePredictionResult(df)


class Distance2SelfResult(AResult):
    """
        Distance2Self prediction result
    """

    def filter_result(self, expressions):
        #TODO: has to be implemented
        pass

    def merge_results(self, others):
        #TODO: has to be implemented
        pass


class CleavageSitePredictionResult(AResult):
    """
        Cleavage site prediction result
    """

    def filter_result(self, expressions):
        """
        Filters a result data frame based on a specified expression consisting
        of a list of triple with (method_name, comparator, threshold). The expression is applied to each row. If any of the
        columns full fill the criteria the row remains.

        :param list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator, threshold)
        :return: CleavageSitePredictionResult - A new filtered result object
        """
        if isinstance(expressions, tuple):
            expressions = [expressions]

        #builde logical expression
        masks = [list(comp(self.loc[:, method], thr)) for method, comp, thr in expressions]

        if len(masks) > 1:
            masks = numpy.logical_and(*masks)
        else:
            masks = masks[0]
        #apply to all rows

        return CleavageSitePredictionResult(self.loc[masks, :])

    def merge_results(self, others):
        """
        Merges results of type CleavageSitePredictionResult and returns the merged result

        :param CleavageSitePredictionResult other: Another CleavageSitePredictionResult object of the same class
        :return: CleavageSitePredictionResult - A new merged CleavageSitePredictionResult object
        """
        if type(others) == type(self):
            others = [others]
        df = self

        for i in xrange(len(others)):
            o = others[i]
            df1a, df2a = df.align(o,)

            o_diff = o.index.difference(df.index)
            d_diff = df.index.difference(o.index)

            if len(d_diff) and len(o_diff):
                df2a.loc[d_diff, "Seq"] = ""
                df1a.loc[o_diff, "Seq"] = ""
            elif len(o_diff):
                df2a.loc[df.index.intersection(o.index), "Seq"] = ""
                df1a.loc[o_diff, "Seq"] = ""
            elif len(d_diff):
                df2a.loc[d_diff, "Seq"] = ""
                df1a.loc[o.index.intersection(df.index), "Seq"] = ""
            else:
                df2a.loc[o.index, "Seq"] = ""

            zero1 = df1a == 0
            zero2 = df2a == 0
            true_zero = zero1 | zero2

            df1 = df1a.fillna(0)
            df2 = df2a.fillna(0)

            df_merged = df1+df2
            false_zero = df_merged == 0
            zero = true_zero & false_zero

            nans = ~true_zero & false_zero
            df_merged = df_merged.where(~zero, other=0)
            df_merged = df_merged.where(~nans, other=numpy.NaN)
            df = df_merged
        return CleavageSitePredictionResult(df)


class CleavageFragmentPredictionResult(AResult):
    """
        Cleavage fragment prediction result
    """

    def filter_result(self, expressions):
        """
        Filters a result data frame based on a specified expression consisting
        of a list of triple with (method_name, comparator, threshold). The expression is applied to each row. If any of the
        columns full fill the criteria the row remains.

        :param list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator, threshold)
        :return: CleavageFragmentPredictionResult - A new filtered result object
        """

        if isinstance(expressions, tuple):
            expressions = [expressions]

        masks = [list(comp(self.loc[:, method], thr)) for method, comp, thr in expressions]

        if len(masks) > 1:
            masks = numpy.logical_and(*masks)
        else:
            masks = masks[0]
        #apply to all rows
        return CleavageFragmentPredictionResult(self.loc[masks, :])

    def merge_results(self, others):
        """
        Merges results of type CleavageFragmentPredictionResult and returns the merged result

        :param CleavageFragmentPredictionResult other: Another CleavageFragmentPredictionResult object of the same class
        :return: CleavageFragmentPredictionResult - A new merged CleavageFragmentPredictionResult object
        """
        if type(others) == type(self):
            others = [others]

        return CleavageFragmentPredictionResult(pandas.concat([self]+others, axis=1))


class TAPPredictionResult(AResult):
    """
        TAP prediction result
    """

    def filter_result(self, expressions):
        """
        filters a result data frame based on a specified expression consisting
        of a list of triple with (method_name, comparator, threshold). The expression is applied to each row. If any of the
        columns full fill the criteria the row remains.

        :param list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator, threshold)
        :return: TAPPredictionResult - A new filtered result object
        """
        if isinstance(expressions, tuple):
            expressions = [expressions]

        masks = [list(comp(self.loc[:, method], thr)) for method, comp, thr in expressions]

        if len(masks) > 1:
            masks = numpy.logical_and(*masks)
        else:
            masks = masks[0]
        #apply to all rows

        return TAPPredictionResult(self.loc[masks, :])

    def merge_results(self, others):
        """
        Merges results of type TAPPredictionResult and returns the merged result

        :param TAPPredictionResult other: Another TAPPredictionResult object of the same class
        :return: TAPPredictionResult - A new merged TAPPredictionResult object
        """
        if type(others) == type(self):
            others = [others]

        return TAPPredictionResult(pandas.concat([self]+others, axis=1))
