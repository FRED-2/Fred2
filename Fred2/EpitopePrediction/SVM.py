# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: EpitopePrediction.SVM
   :synopsis: This module contains all SVM-based epitope prediction tools
.. moduleauthor:: schubert

"""
import svmlight
import collections
import itertools
import pandas
import os

from Fred2.Core.Base import AEpitopePrediction, ASVM
from Fred2.Core.Result import EpitopePredictionResult


class ASVMEpitopePrediction(AEpitopePrediction, ASVM):
    """
        implements default prediction routine for SVM based epitope prediction tools
    """

    def predict(self, peptides, alleles=None, **kwargs):

        if isinstance(peptides, collections.Iterable):
            pep_seqs = {str(p):p for p in peptides}
        else:
            pep_seqs = {str(peptides):peptides}

        if alleles is None:
            allales_string = {conv_a:a for conv_a, a in itertools.izip(self.convert_alleles(self.supportedAlleles),
                                                                       self.supportedAlleles)}
        else:
            allales_string ={conv_a:a.name for conv_a, a in itertools.izip(self.convert_alleles(alleles),alleles)}

        #group peptides by length and
        result = {}
        for length, peps in itertools.groupby(pep_seqs.iterkeys(), key= lambda x: len(x)):
            #load svm model
            encoding = self.encode(peps)

            for a in allales_string.keys():
                model_path = os.path.abspath("../Data/svms/%s/%s_%i"%(self.name, a, length))
                print model_path
                model = svmlight.read_model(model_path)
                pred = svmlight.classify(model, encoding.values())
                result[allales_string[a]] = {}
                for pep, score in itertools.izip(encoding.keys(), pred):
                    result[allales_string[a]][pep_seqs[pep]] = score

        df_result = EpitopePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i, self.name)) for i in df_result.index],
                                                        names=['Seq', 'Method'])
        return df_result


class SVMHC(ASVMEpitopePrediction):
    """
    Implements SVMHC epitope prediction for MHC-I alleles

    """
    __name = "svmhc"
    __alleles = ['HLA-A*02:01', 'HLA-A*02:01', 'HLA-A*11:01', 'HLA-A*11:01', 'HLA-A*24:02', 'HLA-B*15:01',
                 'HLA-B*15:01', 'HLA-B*18:01', 'HLA-B*18:01', 'HLA-B*27:05', 'HLA-B*35:01', 'HLA-B*37:01',
                 'HLA-B*51:01', 'HLA-B*51:01', 'HLA-C*04:01']
    __lengths = [8, 9, 10]

    @property
    def name(self):
        return self.__name

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__lengths

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

    def convert_alleles(self, alleles):
        return ["%s_%s%s"%(a.locus, a.supertype, a.subtype) for a in alleles]

