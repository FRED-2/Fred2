# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: TAPPrediction.PSSM
   :synopsis: This module contains all PSSM-based TAP prediction tools
.. moduleauthor:: schubert

"""

import itertools
import warnings

from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import ATAPPrediction
from Fred2.Core.Result import TAPPredictionResult


class APSSMTAPPrediction(ATAPPrediction):

    def predict(self, peptides, **kwargs):

        def __load_model(length):
            model = "%s_%i"%(self.name, length)
            return getattr( __import__("Fred2.Data.TAPPSSMMatrices", fromlist=[model]), model)


        if isinstance(peptides, Peptide):
            pep_seqs = {str(peptides):peptides}
        else:
            if any(not isinstance(p, Peptide) for p in peptides):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(p):p for p in peptides}

        result = {self.name:{}}
        for length, peps in itertools.groupby(pep_seqs.iterkeys(), key= lambda x: len(x)):
            try:
                pssm = __load_model(length)
            except ImportError:
                    warnings.warn("No model found for %s with length %i"%(self.name, length))
                    continue

            for p in peps:
                score = sum(pssm[i].get(aa, 0.0) for i, aa in enumerate(p))+pssm.get(-1,{}).get("con", 0)
                result[self.name][pep_seqs[p]] = score

        if not result[self.name]:
            raise ValueError("No predictions could be made for given input.")
        df_result = TAPPredictionResult.from_dict(result)

        return df_result


class TAPDoytchinova(APSSMTAPPrediction):
    """
        Implements the TAP prediction model from Doytchinova

        Doytchinova, I., Hemsley, S. and Flower, D. R.
        Transporter associated with antigen processing preselection of peptides binding to the MHC: a bioinformatic evaluation.
        J Immunol, 2004, 173, 6813-6819
    """

    __name = "doytchinova"
    __supported_length = frozenset([9])
    __version = "1.0"

    @property
    def version(self):
        return self.__version

    @property
    def name(self):
        return self.__name

    @property
    def supportedLength(self):
        return self.__supported_length

    def predict(self, peptides, **kwargs):
        return super(TAPDoytchinova, self).predict(peptides, **kwargs)


class SMMTAP(APSSMTAPPrediction):
    """
        Implementation of
        Peters, B., Bulik, S., Tampe, R., Van Endert, P. M., & Holzhuetter, H. G. (2003). Identifying MHC class I
        epitopes by predicting the TAP transport efficiency of epitope precursors. The Journal of Immunology,
        171(4), 1741-1749.
    """
    __name = "smmtap"
    __supported_length = frozenset([9])
    __version = "1.0"

    @property
    def version(self):
        return self.__version

    @property
    def name(self):
        return self.__name

    @property
    def supportedLength(self):
        return self.__supported_length

    def predict(self, peptides, **kwargs):

        def __load_model(length):
            model = "%s_%i"%(self.name, length)
            return getattr(__import__("Fred2.Data.TAPPSSMMatrices", fromlist=[model]), model)


        if isinstance(peptides, Peptide):
            pep_seqs = {str(peptides):peptides}
        else:
            if any(not isinstance(p, Peptide) for p in peptides):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(p):p for p in peptides}

        result = {self.name:{}}
        for length, peps in itertools.groupby(pep_seqs.iterkeys(), key= lambda x: len(x)):
            if length < 9:
                warnings.warn("No model found for %s with length %i"%(self.name, length))
                continue

            try:
                pssm = __load_model(9)
            except ImportError:
                    warnings.warn("No model found for %s with length %i"%(self.name, length))
                    continue

            for p in peps:
                if length <= 9:
                    score = sum(pssm[i].get(aa, 0.0) for i, aa in enumerate(p))
                else:
                    score = sum(pssm[i].get(p[i], 0.0) for i in xrange(3))+pssm[8].get(p[-1], 0.0)
                result[self.name][pep_seqs[p]] = score

        if not result[self.name]:
            raise ValueError("No predictions could be made for given input.")
        df_result = TAPPredictionResult.from_dict(result)

        return df_result