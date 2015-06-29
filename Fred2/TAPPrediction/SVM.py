# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: TAPPrediction.SVM
   :synopsis: This module contains all SVM-based TAP prediction tools
.. moduleauthor:: schubert

"""

import svmlight
import collections
import itertools
import warnings
import pkg_resources

from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import ATAPPrediction, ASVM
from Fred2.Core.Result import TAPPredictionResult


class ASVMTAPPrediction(ATAPPrediction, ASVM):

    def predict(self, peptides,  **kwargs):

        if isinstance(peptides, Peptide):
            pep_seqs = {str(peptides):peptides}
        else:
            if any(not isinstance(p, Peptide) for p in peptides):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(p):p for p in peptides}

        #group peptides by length and

        result = {self.name:{}}
        for length, peps in itertools.groupby(pep_seqs.iterkeys(), key= lambda x: len(x)):
            #load svm model
            encoding = self.encode(peps)

            if length not in self.supportedLength:
                warnings.warn("No model exists for peptides of length %i. Allowed lengths are (%s)"%(length,
                                                                                    ", ".join(self.supportedLength)))
                continue
            model_path = pkg_resources.resource_string("Fred2.Data.svms.%s"%self.name, "%s_%i"%(self.name, length))
            model = svmlight.read_model(model_path)


            pred = svmlight.classify(model, encoding.values())
            result[self.name] = {}
            for pep, score in itertools.izip(encoding.keys(), pred):
                    result[self.name][pep_seqs[pep]] = score

        if not result:
            raise ValueError("No predictions could be made for given input.")
        df_result = TAPPredictionResult.from_dict(result)

        return df_result


class SVMTAP(ASVMTAPPrediction):
    """
        Implements SVMTAP prediction of Doeness et al.

        An SVM method for prediction of TAP affinities.
        Doennes, P. and Kohlbacher, O.
        Integrated modeling of the major events in the MHC class I antigen processing pathway.
        Protein Sci, 2005
    """

    __name = "svmtap"
    __length = [9]

    @property
    def name(self):
        return self.__name

    @property
    def supportedLength(self):
        return self.__length

    def encode(self, peptides):
        """
        Here implements a binary sparse encoding of the peptide

        :param peptides:
        :return: dict(peptide, (tuple(int, list(tuple(int,float)))) -- dictionary with peptide
                 as key and feature encoding as value (see svmlight encoding scheme http://svmlight.joachims.org/)
        """
        AA = {'A': 1, 'C': 2, 'E': 4, 'D': 3, 'G': 6, 'F': 5, 'I': 8, 'H': 7, 'K': 9, 'M': 11, 'L': 10, 'N': 12,
              'Q': 14, 'P': 13, 'S': 16, 'R': 15, 'T': 17, 'W': 19, 'V': 18, 'Y': 20}

        def __encode(pep):
            encoding = []
            offset = 0
            pep_str = str(pep)
            for aa in pep_str:
                    encoding.append((AA[aa]+offset, 1))
                    offset += 20
            return 0, encoding

        if isinstance(peptides, collections.Iterable):
            return {p:__encode(p) for p in peptides}
        else:
            return {peptides:__encode(peptides)}

    def predict(self, peptides,  **kwargs):
       return super(SVMTAP, self).predict(peptides, **kwargs)