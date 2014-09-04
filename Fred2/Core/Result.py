# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Variant
   :synopsis: Contains relevant classes describing results of predictions.
.. moduleauthor:: schubert,

"""
__author__ = 'schubert'

import abc
from pandas import DataFrame


class Result(DataFrame):
    """
        A Result object is a DataFrame with with multi-indexing, where column Id are the predicton model (i.e HLA allele
        for epitope prediction), row ID the target of the prediction (i.e. peptide) and the second row ID the predictor
        (i.e
        the elements for which predictions were run.
    """
