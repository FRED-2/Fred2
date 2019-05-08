# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.AResult
   :synopsis: Contains relevant classes describing results of predictions.
.. moduleauthor:: schubert

"""
__author__ = 'schubert'

import abc
import numpy
import pandas


class AResult(pandas.DataFrame, metaclass=abc.ABCMeta):
    """
        A :class:`~Fred2.Core.Result.AResult` object is a :class:`pandas.DataFrame` with with multi-indexing.

        This class is used as interface and can be extended with custom short-cuts for the sometimes often tedious
        calls in pandas
    """

    @abc.abstractmethod
    def filter_result(self, expressions):
        """
        Filter result based on a list of expressions

        :param list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator,
                                                         threshold)
        :return: A new filtered AResult object
        :rtype: :class:`~Fred2.Core.Result.AResult`
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def merge_results(self, others):
        """
        Merges results of the same type and returns a merged result

        :param others: A (list of) :class:`~Fred2.Core.Result.AResult` object(s) of the same class
        :type others: list(:class:`~Fred2.Core.Result.AResult`)/:class:`~Fred2.Core.Result.AResult`
        :return: A new merged :class:`~Fred2.Core.Result.AResult` object
        :rtype: :class:`~Fred2.Core.Result.AResult`
        """
        raise NotImplementedError()


class EpitopePredictionResult(AResult):
    """
        A :class:`~Fred2.Core.Result.EpitopePredictionResult` object is a DataFrame with multi-indexing, where column
        Ids are the prediction model (i.e HLA :class:`~Fred2.Core.Allele.Allele` for epitope prediction), row ID the
        target of the prediction (i.e. :class:`~Fred2.Core.Peptide.Peptide`) and the second row ID the predictor
        (e.g. BIMAS)

        EpitopePredictionResult

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
        Filters a result data frame based on a specified expression consisting of a list of triple with
        (method_name, comparator, threshold). The expression is applied to each row. If any of the columns fulfill
        the criteria the row remains.

        :param list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator,
                                                        threshold)
        :return: Filtered result object
        :rtype: :class:`~Fred2.Core.Result.EpitopePredictionResult`
        """
        if isinstance(expressions, tuple):
            expressions = [expressions]

        #builde logical expression
        masks = list(map(list, [comp(self.loc[(slice(None), method), :], thr).any(axis=1)
                    for method, comp, thr in expressions]))
        if len(masks) > 1:
            masks = numpy.logical_and(*masks)
        else:
            masks = masks[0]

        #apply to all rows
        idx = [f for f in masks
               for _ in range(len(self.index.levels[1]))]
        return EpitopePredictionResult(self.loc[idx, :])

    def merge_results(self, others):
        """
        Merges results of type :class:`~Fred2.Core.Result.EpitopePredictionResult` and returns the merged result

        :param others: Another (list of) :class:`~Fred2.Core.Result.EpitopePredictionResult`(s)
        :type others: list(:class:`~Fred2.Core.Result.EpitopePredictionResult`)/:class:`~Fred2.Core.Result.EpitopePredictionResult`
        :return: A new merged :class:`~Fred2.Core.Result.EpitopePredictionResult` object
        :rtype: :class:`~Fred2.Core.Result.EpitopePredictionResult`
        """
        df = self.copy(deep=False)

        if type(others) == type(self):
            others = [others]

        for i in range(len(others)):
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
        A :class:`~Fred2.Core.Result.CleavageSitePredictionResult` object is a :class:`pandas.DataFrame` with
        multi-indexing, where column Ids are the prediction scores fo the different prediction methods, as well as the
        amino acid a a specific position, row ID the :class:`~Fred2.Core.Protein.Protein` ID and the position of the
        sequence (starting at 0).



        CleavageSitePredictionResult:

        +--------------+-------------+-------------+-------------+
        | ID           | Pos         |     Seq     | Method_name |
        +==============+=============+=============+=============+
        | protein_ID   |     0       |     S       |     0.56    |
        +              +-------------+-------------+-------------+
        |              |     1       |     Y       |      15     |
        +              +-------------+-------------+-------------+
        |              |     2       |     F       |     0.36    |
        +              +-------------+-------------+-------------+
        |              |     3       |     P       |      10     |
        +--------------+-------------+-------------+-------------+
    """

    def filter_result(self, expressions):
        """
        Filters a result data frame based on a specified expression consisting
        of a list of triple with (method_name, comparator, threshold). The expression is applied to each row. If any of
        the columns fulfill the criteria the row remains.

        :param list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator,
                                                         threshold)
        :return: A new filtered result object
        :rtype: :class:`~Fred2.Core.Result.CleavageSitePredictionResult`
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
        Merges results of type :class:`~Fred2.Core.Result.CleavageSitePredictionResult` and returns the merged result

        :param others: A (list of) :class:`~Fred2.Core.Result.CleavageSitePredictionResult` object(s)
        :type others: list(:class:`~Fred2.Core.Result.CleavageSitePredictionResult`) or
                      :class:`~Fred2.Core.Result.CleavageSitePredictionResult`
        :return: A new merged :class:`~Fred2.Core.Result.CleavageSitePredictionResult` object
        :rtype: :class:`~Fred2.Core.Result.CleavageSitePredictionResult`
        """
        if type(others) == type(self):
            others = [others]
        df = self

        for i in range(len(others)):
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
        A :class:`~Fred2.Core.Result.CleavageFragmentPredictionResult` object is a :class:`pandas.DataFrame` with
        single-indexing, where column  Ids are the prediction scores fo the different prediction methods, and row ID
        the :class:`~Fred2.Core.Peptide.Peptide` object.

        CleavageFragmentPredictionResult:

        +--------------+-------------+
        | Peptide Obj  | Method Name |
        +==============+=============+
        | Peptide1     | -15.34      |
        +--------------+-------------+
        | Peptide2     | 23.34       |
        +--------------+-------------+
    """

    def filter_result(self, expressions):
        """
        Filters a result data frame based on a specified expression consisting of a list of triple with
        (method_name, comparator, threshold). The expression is applied to each row. If any of the columns fulfill the
        criteria the row remains.

        :param list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator,
                                                         threshold)
        :return: A new filtered result object
        :rtype: :class:`~Fred2.Core.Result.CleavageFragmentPredictionResult`
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
        Merges results of type :class:`~Fred2.Core.Result.CleavageFragmentPredictionResult` and returns the merged
        result

        :param others: A (list of) :class:`~Fred2.Core.Result.CleavageFragmentPredictionResult` object(s)
        :type others: list(:class:`~Fred2.Core.Result.CleavageFragmentPredictionResult`) or
                      :class:`~Fred2.Core.Result.CleavageFragmentPredictionResult`
        :return: new merged :class:`~Fred2.Core.Result.CleavageFragmentPredictionResult` object
        :rtype: :class:`~Fred2.Core.Result.CleavageFragmentPredictionResult`
        """
        if type(others) == type(self):
            others = [others]

        return CleavageFragmentPredictionResult(pandas.concat([self]+others, axis=1))


class TAPPredictionResult(AResult):
    """
        A :class:`~Fred2.Core.Result.TAPPredictionResult` object is a :class:`pandas.DataFrame` with single-indexing,
        where column Ids are the ` prediction names of the different prediction methods, and row ID the
        :class:`~Fred2.Core.Peptide.Peptide` object

        TAPPredictionResult:

        +--------------+-------------+
        | Peptide Obj  | Method Name |
        +==============+=============+
        | Peptide1     | -15.34      |
        +--------------+-------------+
        | Peptide2     | 23.34       |
        +--------------+-------------+
    """

    def filter_result(self, expressions):
        """
        Filters a result data frame based on a specified expression consisting
        of a list of triple with (method_name, comparator, threshold). The expression is applied to each row. If any of
        the columns fulfill the criteria the row remains.

        :param list((str,comparator,float)) expressions: A list of triples consisting of (method_name, comparator,
                                                         threshold)
        :return: A new filtered result object
        :rtype: :class:`~Fred2.Core.Result.TAPPredictionResult`
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
        Merges results of type :class:`~Fred2.Core.Result.TAPPredictionResult and returns the merged result

        :param others: A (list of) :class:`~Fred2.Core.Result.TAPPredictionResult` object(s)
        :type others: list(:class:`~Fred2.Core.Result.TAPPredictionResult`) or
                      :class:`~Fred2.Core.Result.TAPPredictionResult`
        :return: A new merged :class:`~Fred2.Core.Result.TAPPredictionResult` object
        :rtype: :class:`~Fred2.Core.Result.TAPPredictionResult``
        """
        if type(others) == type(self):
            others = [others]

        return TAPPredictionResult(pandas.concat([self]+others, axis=1))
