"""
.. module:: TAPPrediction.SVM
   :synopsis: This module contains all SVM-based TAP prediction tools
.. moduleauthor:: schubert

"""
__author__ = 'schubert'

import collections
import itertools
import os
import warnings

from Fred2.Core.Base import ATAPPrediction
from Fred2.Core.Result import TAPPredictionResult


class APSSMTAPPrediction(ATAPPrediction):

    def predict(self, peptides, **kwargs):

        def __load_model(length):
            model = "%s_%i"%(self.name, length)

            #TODO: what if there exists no allele model for this length?
            return getattr( __import__("Fred2.Data.TAPPSSMMatrices.py", fromlist=[model]), model)


        if isinstance(peptides, collections.Iterable):
            pep_seqs = {str(p):p for p in peptides}
        else:
            pep_seqs = {str(peptides):peptides}


        #group peptides by length and

        result = {self.name:{}}
        for length, peps in itertools.groupby(pep_seqs.iterkeys(), key= lambda x: len(x)):
            try:
                pssm = __load_model(length)
            except ImportError:
                    warnings.warn("No model found for %s with length %i"%(self.name, length))
                    continue

            result = {self.name:{}}
            for p in peps:
                score = sum(pssm[i].get(aa, 0.0) for i, aa in enumerate(p))
                result[self.name][pep_seqs[p]] = score

        df_result = TAPPredictionResult.from_dict(result)

        return df_result


class TAPDoytchinova(ATAPPrediction):
    """
        Implements the TAP prediction model from Doytchinova

        Doytchinova, I., Hemsley, S. and Flower, D. R.
        Transporter associated with antigen processing preselection of peptides binding to the MHC: a bioinformatic evaluation.
        J Immunol, 2004, 173, 6813-6819
    """

    __name = "doytchinova"
    __supported_length = [9]

    @property
    def name(self):
        return self.__name

    @property
    def supportedLength(self):
        return self.__supported_length

    def predict(self, peptides, **kwargs):
        return super(TAPDoytchinova, self).predict(peptides, **kwargs)