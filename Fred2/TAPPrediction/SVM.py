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
        """
        Returns TAP predictions for given :class:`~Fred2.Core.Peptide.Peptide`.

        :param peptides: A single :class:`~Fred2.Core.Peptide.Peptide` or a list of :class:`~Fred2.Core.Peptide.Peptide`
        :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`) or :class:`~Fred2.Core.Peptide.Peptide`
        :return: Returns a :class:`~Fred2.Core.Result.TAPPredictionResult` object with the prediction results
        :rtype: :class:`~Fred2.Core.Result.TAPPredictionResult`
        """
        if isinstance(peptides, Peptide):
            pep_seqs = {str(peptides):peptides}
        else:
            pep_seqs = {}
            for p in peptides:
                if not isinstance(p, Peptide):
                    raise ValueError("Input is not of type Protein or Peptide")
                pep_seqs[str(p)] = p

        #group peptides by length and
        chunksize = len(pep_seqs)
        if 'chunks' in kwargs:
            chunksize = kwargs['chunks']

        result = {self.name: {}}
        pep_groups = list(pep_seqs.keys())
        pep_groups.sort(key=len)
        for length, peps in itertools.groupby(pep_groups, key=len):
            #load svm model
            if length not in self.supportedLength:
                warnings.warn("Peptide length of %i is not supported by %s"%(length,self.name))
                continue

            peps = list(peps)
            for i in range(0, len(peps), chunksize):
                encoding = self.encode(peps[i:i+chunksize])

                model_path = pkg_resources.resource_filename("Fred2.Data.svms.%s"%self.name, "%s_%i"%(self.name, length))
                model = svmlight.read_model(model_path)

                pred = svmlight.classify(model, list(encoding.values()))
                for pep, score in zip(list(encoding.keys()), pred):
                        result[self.name][pep_seqs[pep]] = score

        if not result[self.name]:
            raise ValueError("No predictions could be made with "+self.name+" for given input.")
        df_result = TAPPredictionResult.from_dict(result)

        return df_result


class SVMTAP(ASVMTAPPrediction):
    """
    Implements SVMTAP prediction of Doeness et al.

    .. note::

        Doennes, P. and Kohlbacher, O. Integrated modeling of the major events in the MHC class
        I antigen processing pathway. Protein Sci, 2005
    """

    __name = "svmtap"
    __length = frozenset([9])
    __version = "1.0"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    @property
    def supportedLength(self):
        """A list of supported peptide lengths"""
        return self.__length

    def encode(self, peptides):
        """
        Encodes the :class:`~Fred2.Core.Peptide.Peptide` with a binary sparse encoding

        :param list(str) peptides: A list of :class:`~Fred2.Core.Peptide.Peptide`
        :return: Dictionary with :class:`~Fred2.Core.Peptide.Peptide` as key and feature encoding as value (see svmlight
                 encoding scheme http://svmlight.joachims.org/)
        :rtype: dict(:class:`~Fred2.Core.Peptide.Peptide`, (tuple(int, list(tuple(int,float))))
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
            return {p: __encode(p) for p in peptides}
        else:
            return {peptides: __encode(peptides)}

    def predict(self, peptides, **kwargs):
        """
        Returns predictions for given :class:`~Fred2.Core.Peptide.Peptide`.

        :param peptides: A single :class:`~Fred2.Core.Peptide.Peptide` or a list of :class:`~Fred2.Core.Peptide.Peptide`
        :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`) or :class:`~Fred2.Core.Peptide.Peptide`
        :return: Returns a :class:`~Fred2.Core.Result.TAPPredictionResult` object with the prediction results
        :rtype: :class:`~Fred2.Core.Result.TAPPredictionResult`
        """
        return super(SVMTAP, self).predict(peptides, **kwargs)