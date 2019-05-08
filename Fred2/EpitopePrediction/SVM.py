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
import os
import warnings

import pandas
import pkg_resources

from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import AEpitopePrediction, ASVM
from Fred2.Core.Result import EpitopePredictionResult
from Fred2.Data.svms.unitope.UniTope_encodedAlleles import UniTope_encodedAlleles


class ASVMEpitopePrediction(AEpitopePrediction, ASVM):
    """
        Implements default prediction routine for SVM based epitope prediction tools
    """

    def predict(self, peptides, alleles=None, **kwargs):
        """
        Returns predictions for given peptides an alleles. If no alleles are given, predictions for all available models
        are made.

        :param peptides: A single :class:`~Fred2.Core.Peptide.Peptide` or a list of :class:`~Fred2.Core.Peptide.Peptide`
        :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`) or :class:`~Fred2.Core.Peptide.Peptide`
        :param alleles: A list of :class:`~Fred2.Core.Allele.Allele`
        :type alleles: list(:class:`~Fred2.Core.Allele.Allele`) or :class:`~Fred2.Core.Allele.Allele`
        :param kwargs: optional parameter (not used yet)
        :return: Returns a :class:`~Fred2.Core.Result.EpitopePredictionResult` object with the prediction results
        :rtype: :class:`~Fred2.Core.Result.EpitopePredictionResult`
        """
        if isinstance(peptides, Peptide):
            pep_seqs = {str(peptides): peptides}
        else:
            pep_seqs = {}
            for p in peptides:
                if not isinstance(p, Peptide):
                    raise ValueError("Input is not of type Protein or Peptide")
                pep_seqs[str(p)] = p

        chunksize = len(pep_seqs)
        if 'chunks' in kwargs:
            chunksize = kwargs['chunks']

        if alleles is None:
            al = [Allele(a) for a in self.supportedAlleles]
            allales_string = {conv_a: a for conv_a, a in zip(self.convert_alleles(al), al)}
        else:
            if isinstance(alleles, Allele):
                alleles = [alleles]
            if any(not isinstance(p, Allele) for p in alleles):
                raise ValueError("Input is not of type Allele")
            allales_string = {conv_a: a for conv_a, a in zip(self.convert_alleles(alleles), alleles)}

        # group peptides by length and
        result = {}
        pep_groups = list(pep_seqs.keys())
        pep_groups.sort(key=len)
        for length, peps in itertools.groupby(pep_groups, key=len):
            # load svm model

            if length not in self.supportedLength:
                warnings.warn("Peptide length of %i is not supported by %s" % (length, self.name))
                continue

            peps = list(peps)
            for i in range(0, len(peps), chunksize):
                encoding = self.encode(peps[i:i+chunksize])

                for a in list(allales_string.keys()):
                    model_path = pkg_resources.resource_filename("Fred2.Data.svms.%s" % self.name, "%s_%i" % (a, length))
                    if not os.path.exists(model_path):
                        warnings.warn("No model exists for peptides of length %i or allele %s." % (length,
                                                                                                   allales_string[a].name))
                        continue

                    model = svmlight.read_model(model_path)
                    pred = svmlight.classify(model, list(encoding.values()))
                    if allales_string[a] not in result:
                        result[allales_string[a]] = {}
                    for pep, score in zip(list(encoding.keys()), pred):
                        result[allales_string[a]][pep_seqs[pep]] = score

        if not result:
            raise ValueError("No predictions could be made for given input. Check your "
                             "epitope length and HLA allele combination.")
        df_result = EpitopePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i, self.name)) for i in df_result.index],
                                                        names=['Seq', 'Method'])
        return df_result


class SVMHC(ASVMEpitopePrediction):
    """
    Implements SVMHC epitope prediction for MHC-I alleles (SYFPEITHI models).

    .. note::

        Doennes, P. and Kohlbacher, O. SVMHC: a server for prediction of MHC-binding peptides.
        Nucleic Acids Res, 2006, 34, W194-W197
    """
    __name = "svmhc"
    __alleles = frozenset(['HLA-A*02:01', 'HLA-A*02:01', 'HLA-A*11:01', 'HLA-A*11:01', 'HLA-A*24:02', 'HLA-B*15:01',
                           'HLA-B*15:01', 'HLA-B*18:01', 'HLA-B*18:01', 'HLA-B*27:05', 'HLA-B*35:01', 'HLA-B*37:01',
                           'HLA-B*51:01', 'HLA-B*51:01', 'HLA-C*04:01'])
    __supported_length = frozenset([8, 9, 10])
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
    def supportedAlleles(self):
        """
        A list of supported allele models
        """
        return self.__alleles

    @property
    def supportedLength(self):
        """
        A list of supported :class:`~Fred2.Core.Peptide.Peptide` lengths
        """
        return self.__supported_length

    def encode(self, peptides):
        """
        Encodes the input with binary sparse encoding of the :class:`~Fred2.Core.Peptide.Peptide`

        :param str peptides: A list of :class:`~Fred2.Core.Peptide.Peptide` sequences
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
                encoding.append((AA[aa] + offset, 1))
                offset += 20
            return 0, encoding

        if isinstance(peptides, collections.Iterable):
            return {p: __encode(p) for p in peptides}
        else:
            return {peptides: __encode(peptides)}

    def convert_alleles(self, alleles):
        """
        Converts :class:`~Fred2.Core.Allele.Allele` into the internal :class:`~Fred2.Core.Allele.Allele`
        representation of the predictor and returns a string representation

        :param alleles: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is needed
        :type alleles: list(:class:`~Fred2.Core.Allele.Allele`)
        :return: Returns a string representation of the input :class:`~Fred2.Core.Allele.Allele`
        :rtype: list(str)
        """
        return ["%s_%s%s" % (a.locus, a.supertype, a.subtype) for a in alleles]


class UniTope(ASVMEpitopePrediction):
    """
    Implements UniTope prediction for MHC-I.

    .. note::

        Toussaint, N. C., Feldhahn, M., Ziehm, M., Stevanovic, S., & Kohlbacher, O. (2011, August).
        T-cell epitope prediction based on self-tolerance. In Proceedings of the 2nd ACM Conference on
        Bioinformatics, Computational Biology and Biomedicine (pp. 584-588). ACM.
    """

    __name = "unitope"
    __supported_length = frozenset([8, 9, 10])
    __alleles = frozenset(['HLA-B*51:43', 'HLA-A*11:02', 'HLA-C*07:27', 'HLA-B*35:55', 'HLA-B*49:02', 'HLA-A*11:17',
                           'HLA-B*51:02', 'HLA-B*15:61', 'HLA-A*33:04', 'HLA-C*02:12', 'HLA-B*55:07', 'HLA-A*24:13',
                           'HLA-C*01:02', 'HLA-A*24:15', 'HLA-A*24:59', 'HLA-A*11:19', 'HLA-B*15:63', 'HLA-A*32:06',
                           'HLA-B*55:05', 'HLA-B*47:04', 'HLA-B*67:01', 'HLA-C*16:09', 'HLA-B*53:02', 'HLA-B*56:04',
                           'HLA-B*37:12', 'HLA-B*55:03', 'HLA-B*39:30', 'HLA-C*04:20', 'HLA-B*08:21', 'HLA-A*74:02',
                           'HLA-B*53:04', 'HLA-B*56:06', 'HLA-A*01:06', 'HLA-C*03:29', 'HLA-B*55:01', 'HLA-A*36:03',
                           'HLA-A*02:87', 'HLA-A*24:17', 'HLA-A*02:02', 'HLA-A*24:50', 'HLA-C*03:22', 'HLA-C*12:09',
                           'HLA-A*02:17', 'HLA-B*55:08', 'HLA-B*44:15', 'HLA-B*15:46', 'HLA-C*16:02', 'HLA-C*03:24',
                           'HLA-C*07:44', 'HLA-A*02:19', 'HLA-A*66:01', 'HLA-A*01:08', 'HLA-B*15:44', 'HLA-A*02:77',
                           'HLA-C*18:02', 'HLA-C*07:29', 'HLA-B*44:06', 'HLA-C*07:21', 'HLA-A*66:03', 'HLA-A*74:04',
                           'HLA-A*02:80', 'HLA-B*44:11', 'HLA-B*15:42', 'HLA-A*02:75', 'HLA-B*15:65', 'HLA-B*44:08',
                           'HLA-A*02:97', 'HLA-A*24:52', 'HLA-A*66:05', 'HLA-A*74:06', 'HLA-C*07:25', 'HLA-B*27:02',
                           'HLA-B*15:40', 'HLA-A*02:73', 'HLA-B*40:05', 'HLA-B*15:67', 'HLA-B*58:05', 'HLA-A*33:06',
                           'HLA-B*08:03', 'HLA-C*03:26', 'HLA-B*15:69', 'HLA-A*30:08', 'HLA-C*07:09', 'HLA-C*06:11',
                           'HLA-B*18:21', 'HLA-B*41:07', 'HLA-B*15:91', 'HLA-C*15:10', 'HLA-C*07:34', 'HLA-C*12:06',
                           'HLA-B*35:39', 'HLA-C*06:13', 'HLA-B*08:27', 'HLA-B*38:01', 'HLA-A*31:02', 'HLA-C*15:12',
                           'HLA-B*27:08', 'HLA-A*74:08', 'HLA-B*13:15', 'HLA-A*02:79', 'HLA-B*53:06', 'HLA-A*25:05',
                           'HLA-A*03:18', 'HLA-B*95:08', 'HLA-B*15:83', 'HLA-A*24:54', 'HLA-B*53:11', 'HLA-C*08:11',
                           'HLA-B*15:05', 'HLA-B*53:08', 'HLA-B*56:02', 'HLA-A*24:66', 'HLA-B*38:05', 'HLA-B*15:81',
                           'HLA-A*24:56', 'HLA-A*33:09', 'HLA-B*51:38', 'HLA-B*15:01', 'HLA-B*35:58', 'HLA-B*15:48',
                           'HLA-B*35:17', 'HLA-B*54:03', 'HLA-C*03:04', 'HLA-A*02:71', 'HLA-B*55:22', 'HLA-B*51:14',
                           'HLA-B*95:06', 'HLA-A*68:29', 'HLA-B*48:06', 'HLA-B*95:21', 'HLA-B*18:01', 'HLA-C*03:02',
                           'HLA-A*02:36', 'HLA-B*55:09', 'HLA-B*51:16', 'HLA-B*95:04', 'HLA-A*24:58', 'HLA-A*03:14',
                           'HLA-A*02:39', 'HLA-B*35:25', 'HLA-C*12:04', 'HLA-B*40:25', 'HLA-B*51:18', 'HLA-B*35:50',
                           'HLA-C*07:22', 'HLA-B*15:10', 'HLA-A*03:16', 'HLA-A*26:14', 'HLA-A*31:04', 'HLA-C*15:14',
                           'HLA-B*35:15', 'HLA-B*44:48', 'HLA-B*27:28', 'HLA-B*39:36', 'HLA-A*03:10', 'HLA-A*31:06',
                           'HLA-C*15:16', 'HLA-A*02:21', 'HLA-A*33:01', 'HLA-B*15:87', 'HLA-B*40:53', 'HLA-B*35:36',
                           'HLA-B*18:07', 'HLA-A*02:42', 'HLA-B*44:40', 'HLA-B*15:12', 'HLA-A*30:06', 'HLA-B*27:05',
                           'HLA-B*35:30', 'HLA-C*07:06', 'HLA-C*05:11', 'HLA-B*78:02', 'HLA-B*44:42', 'HLA-B*07:37',
                           'HLA-B*27:07', 'HLA-B*40:21', 'HLA-B*38:07', 'HLA-A*01:09', 'HLA-B*18:03', 'HLA-C*07:24',
                           'HLA-A*02:29', 'HLA-A*30:02', 'HLA-B*27:01', 'HLA-C*07:45', 'HLA-A*26:04', 'HLA-B*95:02',
                           'HLA-B*15:85', 'HLA-B*40:51', 'HLA-A*74:01', 'HLA-B*35:11', 'HLA-C*03:14', 'HLA-C*03:06',
                           'HLA-B*27:03', 'HLA-C*07:43', 'HLA-B*55:20', 'HLA-A*26:06', 'HLA-B*38:14', 'HLA-B*08:22',
                           'HLA-B*35:32', 'HLA-C*07:04', 'HLA-A*26:08', 'HLA-B*08:20', 'HLA-C*16:04', 'HLA-A*03:12',
                           'HLA-B*51:32', 'HLA-B*35:19', 'HLA-B*07:02', 'HLA-B*18:08', 'HLA-A*24:25', 'HLA-C*16:06',
                           'HLA-A*02:55', 'HLA-A*02:10', 'HLA-B*51:30', 'HLA-A*02:25', 'HLA-A*23:14', 'HLA-A*68:26',
                           'HLA-B*13:10', 'HLA-A*29:11', 'HLA-B*27:18', 'HLA-A*02:57', 'HLA-A*02:12', 'HLA-A*26:22',
                           'HLA-B*40:57', 'HLA-A*02:27', 'HLA-A*23:12', 'HLA-A*02:99', 'HLA-B*44:44', 'HLA-A*02:09',
                           'HLA-A*02:51', 'HLA-A*02:14', 'HLA-A*26:20', 'HLA-B*08:28', 'HLA-B*35:34', 'HLA-C*07:02',
                           'HLA-A*23:10', 'HLA-B*27:12', 'HLA-A*92:09', 'HLA-B*44:46', 'HLA-B*07:33', 'HLA-A*02:07',
                           'HLA-C*17:02', 'HLA-B*58:12', 'HLA-C*02:02', 'HLA-B*40:39', 'HLA-A*33:05', 'HLA-A*74:03',
                           'HLA-C*05:14', 'HLA-B*48:08', 'HLA-B*13:01', 'HLA-B*27:35', 'HLA-B*15:34', 'HLA-B*45:03',
                           'HLA-B*44:18', 'HLA-B*15:19', 'HLA-A*33:03', 'HLA-B*39:22', 'HLA-B*13:03', 'HLA-B*08:26',
                           'HLA-B*15:32', 'HLA-B*45:05', 'HLA-B*07:38', 'HLA-B*53:01', 'HLA-A*02:22', 'HLA-A*92:03',
                           'HLA-B*39:20', 'HLA-C*01:10', 'HLA-B*08:24', 'HLA-C*16:08', 'HLA-B*53:03', 'HLA-B*51:36',
                           'HLA-B*56:15', 'HLA-A*74:09', 'HLA-B*39:26', 'HLA-A*29:13', 'HLA-A*24:21', 'HLA-A*02:69',
                           'HLA-B*53:05', 'HLA-A*26:28', 'HLA-C*03:28', 'HLA-B*15:11', 'HLA-A*26:26', 'HLA-B*35:65Q',
                           'HLA-A*24:05', 'HLA-B*48:07', 'HLA-B*13:02', 'HLA-A*11:25', 'HLA-A*02:05', 'HLA-C*03:21',
                           'HLA-A*02:90', 'HLA-A*24:03', 'HLA-A*02:37', 'HLA-B*50:01', 'HLA-A*02:76', 'HLA-A*02:03',
                           'HLA-C*03:23', 'HLA-B*08:29', 'HLA-A*02:16', 'HLA-A*92:07', 'HLA-A*02:74', 'HLA-C*02:06',
                           'HLA-B*07:26', 'HLA-A*02:01', 'HLA-C*03:25', 'HLA-C*05:10', 'HLA-B*45:07', 'HLA-A*74:05',
                           'HLA-A*02:18', 'HLA-A*26:24', 'HLA-B*08:17', 'HLA-A*92:01', 'HLA-A*02:72', 'HLA-C*02:04',
                           'HLA-A*33:07', 'HLA-C*03:27', 'HLA-B*47:02', 'HLA-A*74:07', 'HLA-B*47:01', 'HLA-B*15:36',
                           'HLA-B*15:17', 'HLA-B*13:09', 'HLA-A*01:17', 'HLA-C*04:16', 'HLA-B*53:07', 'HLA-B*15:04',
                           'HLA-A*30:09', 'HLA-B*40:49', 'HLA-B*40:37', 'HLA-B*18:20', 'HLA-B*07:10', 'HLA-B*14:01',
                           'HLA-B*27:09', 'HLA-C*04:14', 'HLA-B*53:09', 'HLA-B*40:04', 'HLA-A*02:78', 'HLA-B*15:58',
                           'HLA-B*40:31', 'HLA-B*18:22', 'HLA-B*07:12', 'HLA-B*39:24', 'HLA-B*13:06', 'HLA-A*02:92',
                           'HLA-B*39:10', 'HLA-A*03:01', 'HLA-A*24:26', 'HLA-B*40:33', 'HLA-B*18:24', 'HLA-A*32:02',
                           'HLA-B*95:09', 'HLA-B*53:10', 'HLA-B*15:99', 'HLA-A*30:07', 'HLA-B*13:11', 'HLA-A*24:28',
                           'HLA-A*24:07', 'HLA-B*15:39', 'HLA-C*16:01', 'HLA-B*07:28', 'HLA-B*39:28', 'HLA-B*78:05',
                           'HLA-A*33:08', 'HLA-B*52:07', 'HLA-A*11:23', 'HLA-A*26:01', 'HLA-B*35:60', 'HLA-B*95:20',
                           'HLA-B*07:51', 'HLA-B*40:08', 'HLA-C*03:07', 'HLA-B*35:66', 'HLA-B*40:44', 'HLA-B*95:22',
                           'HLA-A*03:15', 'HLA-B*15:74', 'HLA-B*39:06', 'HLA-B*40:28', 'HLA-C*03:05', 'HLA-B*37:07',
                           'HLA-B*51:34', 'HLA-B*95:05', 'HLA-B*13:04', 'HLA-A*26:02', 'HLA-A*03:17', 'HLA-A*01:13',
                           'HLA-C*03:03', 'HLA-C*08:08', 'HLA-B*37:05', 'HLA-A*26:03', 'HLA-A*02:58', 'HLA-B*15:78',
                           'HLA-B*44:28', 'HLA-C*17:04', 'HLA-B*58:14', 'HLA-B*55:12', 'HLA-C*04:18', 'HLA-B*51:19',
                           'HLA-A*24:08', 'HLA-B*15:31', 'HLA-A*30:01', 'HLA-A*29:15', 'HLA-B*27:06', 'HLA-B*15:13',
                           'HLA-C*08:05', 'HLA-C*04:11', 'HLA-B*40:52', 'HLA-A*03:19', 'HLA-B*44:20', 'HLA-B*18:06',
                           'HLA-B*15:55', 'HLA-B*07:36', 'HLA-C*08:07', 'HLA-B*38:08', 'HLA-B*40:50', 'HLA-A*02:52',
                           'HLA-B*78:01', 'HLA-C*14:05', 'HLA-B*15:57', 'HLA-C*05:12', 'HLA-C*06:08', 'HLA-B*37:09',
                           'HLA-B*40:02', 'HLA-B*56:05', 'HLA-B*95:01', 'HLA-B*08:23', 'HLA-A*26:05', 'HLA-B*78:03',
                           'HLA-B*18:02', 'HLA-B*15:51', 'HLA-B*15:72', 'HLA-A*30:03', 'HLA-B*95:07', 'HLA-B*45:02',
                           'HLA-A*26:07', 'HLA-B*35:62', 'HLA-B*15:53', 'HLA-B*15:70', 'HLA-B*39:02', 'HLA-B*48:13',
                           'HLA-C*03:09', 'HLA-B*35:57', 'HLA-B*51:33', 'HLA-A*02:93', 'HLA-A*31:05', 'HLA-B*54:10',
                           'HLA-B*37:01', 'HLA-B*18:09', 'HLA-A*24:22', 'HLA-A*02:96', 'HLA-C*01:11', 'HLA-A*02:95',
                           'HLA-A*80:01', 'HLA-B*15:84', 'HLA-A*24:24', 'HLA-C*02:15', 'HLA-B*15:38', 'HLA-C*12:12',
                           'HLA-A*03:13', 'HLA-A*23:13', 'HLA-B*56:17', 'HLA-A*24:43', 'HLA-B*58:13', 'HLA-C*17:03',
                           'HLA-B*39:41', 'HLA-C*08:09', 'HLA-C*04:15', 'HLA-B*40:56', 'HLA-A*02:54', 'HLA-B*51:31',
                           'HLA-A*31:09', 'HLA-A*92:08', 'HLA-B*15:92', 'HLA-A*24:41', 'HLA-B*27:04', 'HLA-B*07:32',
                           'HLA-C*04:13', 'HLA-B*55:19', 'HLA-B*40:54', 'HLA-A*02:56', 'HLA-A*26:23', 'HLA-A*24:62',
                           'HLA-B*45:04', 'HLA-B*07:39', 'HLA-B*15:35', 'HLA-C*04:03', 'HLA-C*08:02', 'HLA-B*40:38',
                           'HLA-B*52:08', 'HLA-B*48:11', 'HLA-C*05:15', 'HLA-B*15:33', 'HLA-C*12:18', 'HLA-B*40:58',
                           'HLA-A*92:02', 'HLA-B*39:23', 'HLA-A*11:04', 'HLA-C*01:13', 'HLA-B*40:34', 'HLA-B*51:37',
                           'HLA-B*48:15', 'HLA-B*40:12', 'HLA-B*41:01', 'HLA-A*68:28', 'HLA-A*11:22', 'HLA-B*08:25',
                           'HLA-C*07:18', 'HLA-B*40:36', 'HLA-A*26:29', 'HLA-C*07:36', 'HLA-B*35:48', 'HLA-B*39:27',
                           'HLA-A*02:91', 'HLA-A*24:20', 'HLA-A*68:24', 'HLA-B*35:43', 'HLA-A*31:13', 'HLA-C*15:03',
                           'HLA-A*24:04', 'HLA-A*26:21', 'HLA-B*51:03', 'HLA-A*24:29', 'HLA-C*12:15', 'HLA-B*35:45',
                           'HLA-A*11:11', 'HLA-A*30:04', 'HLA-A*24:55', 'HLA-A*31:11', 'HLA-A*24:02', 'HLA-A*26:27',
                           'HLA-A*24:57', 'HLA-B*50:04', 'HLA-B*51:28', 'HLA-A*68:06', 'HLA-B*15:16', 'HLA-B*51:39',
                           'HLA-B*15:03', 'HLA-A*68:13', 'HLA-A*02:28', 'HLA-A*92:04', 'HLA-B*15:96', 'HLA-B*15:82',
                           'HLA-A*68:04', 'HLA-A*32:13', 'HLA-B*15:18', 'HLA-B*15:37', 'HLA-A*68:15', 'HLA-A*92:06',
                           'HLA-B*39:14', 'HLA-B*35:47', 'HLA-C*07:39', 'HLA-B*40:30', 'HLA-B*52:04', 'HLA-B*07:13',
                           'HLA-A*34:01', 'HLA-B*82:02', 'HLA-C*04:19', 'HLA-B*55:13', 'HLA-A*32:03', 'HLA-B*58:11',
                           'HLA-A*34:03', 'HLA-B*50:02', 'HLA-B*51:13', 'HLA-A*01:14', 'HLA-B*13:08', 'HLA-B*15:15',
                           'HLA-C*04:17', 'HLA-A*01:20', 'HLA-A*24:38', 'HLA-C*12:17', 'HLA-B*39:29', 'HLA-B*78:04',
                           'HLA-C*07:31', 'HLA-B*40:07', 'HLA-B*39:38Q', 'HLA-C*15:07', 'HLA-B*40:32', 'HLA-B*35:02',
                           'HLA-C*06:02', 'HLA-B*44:21', 'HLA-B*35:41', 'HLA-B*46:05', 'HLA-B*40:01', 'HLA-B*46:01',
                           'HLA-B*44:51', 'HLA-A*24:06', 'HLA-B*07:17', 'HLA-A*24:27', 'HLA-A*68:02', 'HLA-B*07:50',
                           'HLA-B*40:66', 'HLA-B*38:02', 'HLA-B*55:15', 'HLA-A*68:17', 'HLA-C*07:20', 'HLA-B*13:16',
                           'HLA-B*35:61', 'HLA-B*35:23', 'HLA-C*07:11', 'HLA-B*18:18', 'HLA-B*38:03', 'HLA-A*68:19',
                           'HLA-B*37:06', 'HLA-B*13:14', 'HLA-B*35:67', 'HLA-A*68:33', 'HLA-B*14:02', 'HLA-C*07:40',
                           'HLA-B*57:05', 'HLA-B*51:46', 'HLA-A*29:02', 'HLA-B*55:11', 'HLA-B*54:09', 'HLA-B*37:04',
                           'HLA-B*35:04', 'HLA-B*15:30', 'HLA-B*07:19', 'HLA-B*57:07', 'HLA-A*01:12', 'HLA-B*15:75',
                           'HLA-B*37:02', 'HLA-A*02:38', 'HLA-A*32:07', 'HLA-A*32:12', 'HLA-A*24:64', 'HLA-B*57:01',
                           'HLA-B*44:29', 'HLA-B*07:34', 'HLA-A*24:47', 'HLA-B*35:21', 'HLA-C*07:13', 'HLA-B*15:98',
                           'HLA-B*53:12', 'HLA-B*15:56', 'HLA-C*14:04', 'HLA-A*68:08', 'HLA-B*40:03', 'HLA-A*34:06',
                           'HLA-A*29:14', 'HLA-A*02:30', 'HLA-C*08:06', 'HLA-B*44:24', 'HLA-C*04:10', 'HLA-B*15:50',
                           'HLA-C*14:02', 'HLA-B*52:02', 'HLA-B*39:01', 'HLA-A*02:65', 'HLA-A*24:44', 'HLA-B*37:08',
                           'HLA-B*35:08', 'HLA-C*08:12', 'HLA-C*08:04', 'HLA-B*39:03', 'HLA-A*02:67', 'HLA-B*48:12',
                           'HLA-B*15:73', 'HLA-B*35:29', 'HLA-A*02:34', 'HLA-A*24:23', 'HLA-B*13:12', 'HLA-A*68:31',
                           'HLA-B*44:25', 'HLA-B*40:09', 'HLA-A*02:61', 'HLA-B*48:10', 'HLA-B*15:71', 'HLA-B*35:27',
                           'HLA-C*07:15', 'HLA-A*74:10', 'HLA-B*39:37', 'HLA-B*35:63', 'HLA-A*02:84', 'HLA-B*15:52',
                           'HLA-B*39:34', 'HLA-B*58:01', 'HLA-C*12:13', 'HLA-B*57:03', 'HLA-C*03:15', 'HLA-C*07:35',
                           'HLA-A*24:65', 'HLA-C*06:10', 'HLA-B*58:07', 'HLA-C*12:11', 'HLA-B*40:13', 'HLA-C*03:17',
                           'HLA-B*35:69', 'HLA-A*02:20', 'HLA-B*45:06', 'HLA-A*24:67', 'HLA-B*95:18', 'HLA-A*30:18',
                           'HLA-B*48:02', 'HLA-B*55:18', 'HLA-C*04:12', 'HLA-A*30:12', 'HLA-C*03:19', 'HLA-C*14:08',
                           'HLA-A*24:61', 'HLA-A*66:02', 'HLA-A*34:02', 'HLA-A*29:10', 'HLA-B*40:48', 'HLA-B*56:18',
                           'HLA-B*57:09', 'HLA-A*30:10', 'HLA-A*68:27', 'HLA-B*15:54', 'HLA-C*14:06', 'HLA-A*24:63',
                           'HLA-A*34:04', 'HLA-A*29:12', 'HLA-B*52:01', 'HLA-B*39:08', 'HLA-B*48:16', 'HLA-A*30:16',
                           'HLA-A*11:05', 'HLA-B*52:06', 'HLA-A*31:12', 'HLA-B*07:35', 'HLA-B*51:29', 'HLA-B*56:11',
                           'HLA-C*03:32', 'HLA-C*08:03', 'HLA-C*12:05', 'HLA-C*06:09', 'HLA-A*68:20', 'HLA-B*51:24',
                           'HLA-C*03:30', 'HLA-C*12:19', 'HLA-B*95:16', 'HLA-A*02:63', 'HLA-A*30:13', 'HLA-B*54:07',
                           'HLA-A*02:49', 'HLA-B*52:09', 'HLA-B*40:15', 'HLA-A*11:07', 'HLA-B*07:07', 'HLA-B*51:07',
                           'HLA-B*46:04', 'HLA-B*35:18', 'HLA-B*48:14', 'HLA-B*27:26', 'HLA-A*68:25', 'HLA-A*31:07',
                           'HLA-C*03:34', 'HLA-C*01:12', 'HLA-A*29:16', 'HLA-C*03:13', 'HLA-B*35:49', 'HLA-C*07:37',
                           'HLA-B*27:21', 'HLA-B*35:64', 'HLA-C*15:02', 'HLA-B*46:06', 'HLA-B*58:08', 'HLA-B*07:23',
                           'HLA-C*07:26', 'HLA-B*44:37', 'HLA-C*12:14', 'HLA-B*41:08', 'HLA-A*02:47', 'HLA-B*56:14',
                           'HLA-B*42:07', 'HLA-C*06:03', 'HLA-A*69:01', 'HLA-B*07:25', 'HLA-B*27:25', 'HLA-B*08:11',
                           'HLA-B*95:10', 'HLA-A*68:10', 'HLA-B*40:65', 'HLA-A*32:11Q', 'HLA-B*15:95', 'HLA-B*42:05',
                           'HLA-B*35:42', 'HLA-A*25:01', 'HLA-A*11:08', 'HLA-C*02:03', 'HLA-C*08:01', 'HLA-B*95:12',
                           'HLA-C*02:05', 'HLA-B*41:04', 'HLA-C*02:08', 'HLA-B*56:10', 'HLA-B*40:11', 'HLA-B*35:44',
                           'HLA-B*38:11', 'HLA-A*11:27', 'HLA-B*51:26', 'HLA-B*14:03', 'HLA-B*95:14', 'HLA-A*68:14',
                           'HLA-B*41:06', 'HLA-B*51:23', 'HLA-C*07:38', 'HLA-C*06:05', 'HLA-C*16:07', 'HLA-A*26:12',
                           'HLA-B*27:16', 'HLA-B*18:11', 'HLA-B*27:24', 'HLA-B*42:08', 'HLA-B*35:01', 'HLA-B*07:16',
                           'HLA-A*02:45', 'HLA-B*54:06', 'HLA-A*24:39', 'HLA-B*38:06', 'HLA-A*03:24', 'HLA-B*18:13',
                           'HLA-A*03:06', 'HLA-B*35:03', 'HLA-C*07:19', 'HLA-B*82:01', 'HLA-A*11:09', 'HLA-A*29:04',
                           'HLA-B*81:01', 'HLA-B*44:50', 'HLA-A*03:04', 'HLA-C*15:06', 'HLA-B*46:02', 'HLA-B*58:04',
                           'HLA-A*11:14', 'HLA-B*38:13', 'HLA-C*12:10', 'HLA-B*08:14', 'HLA-B*40:63', 'HLA-B*27:32',
                           'HLA-A*03:02', 'HLA-C*15:04', 'HLA-B*07:21', 'HLA-C*12:16', 'HLA-B*40:61', 'HLA-A*02:41',
                           'HLA-C*07:30', 'HLA-A*32:04', 'HLA-B*13:17', 'HLA-B*27:15', 'HLA-B*52:10', 'HLA-A*34:08',
                           'HLA-B*38:15', 'HLA-C*02:10', 'HLA-B*40:59', 'HLA-B*46:08', 'HLA-C*07:14', 'HLA-B*35:26',
                           'HLA-C*05:06', 'HLA-B*27:13', 'HLA-B*07:03', 'HLA-A*68:12', 'HLA-A*36:02', 'HLA-B*57:04',
                           'HLA-B*08:18', 'HLA-A*03:08', 'HLA-A*68:03', 'HLA-C*07:16', 'HLA-B*35:24', 'HLA-B*08:01',
                           'HLA-A*68:16', 'HLA-B*07:15', 'HLA-B*57:06', 'HLA-B*08:16', 'HLA-B*27:34', 'HLA-C*15:08',
                           'HLA-A*68:01', 'HLA-C*07:10', 'HLA-B*35:22', 'HLA-B*39:32', 'HLA-B*15:97', 'HLA-B*35:05',
                           'HLA-B*81:02', 'HLA-B*41:02', 'HLA-B*44:22', 'HLA-C*03:11', 'HLA-B*40:43', 'HLA-B*27:36',
                           'HLA-C*05:01', 'HLA-C*07:12', 'HLA-B*35:20', 'HLA-B*08:05', 'HLA-B*35:07', 'HLA-B*41:05',
                           'HLA-B*44:12', 'HLA-B*27:11', 'HLA-B*35:09', 'HLA-B*14:06', 'HLA-C*12:08', 'HLA-B*07:44',
                           'HLA-B*14:05', 'HLA-B*27:27', 'HLA-B*57:08', 'HLA-A*24:31', 'HLA-B*51:08', 'HLA-B*15:80',
                           'HLA-B*59:01', 'HLA-A*34:07', 'HLA-B*55:14', 'HLA-B*44:35', 'HLA-B*15:86', 'HLA-A*24:37',
                           'HLA-C*04:05', 'HLA-B*67:02', 'HLA-A*31:10', 'HLA-A*02:64', 'HLA-B*56:13', 'HLA-C*05:04',
                           'HLA-C*14:03', 'HLA-B*35:37', 'HLA-B*13:13', 'HLA-B*27:19', 'HLA-A*26:18', 'HLA-B*39:09',
                           'HLA-A*02:66', 'HLA-B*49:04', 'HLA-C*02:14', 'HLA-A*02:35', 'HLA-A*68:09', 'HLA-B*15:76',
                           'HLA-A*74:11', 'HLA-B*27:17', 'HLA-B*15:02', 'HLA-B*39:04', 'HLA-A*02:60', 'HLA-C*02:16',
                           'HLA-B*51:20', 'HLA-B*35:28', 'HLA-B*08:12', 'HLA-B*15:09', 'HLA-B*42:01', 'HLA-B*08:07',
                           'HLA-A*03:20', 'HLA-B*35:46', 'HLA-A*25:04', 'HLA-A*32:10', 'HLA-A*02:59', 'HLA-B*35:68',
                           'HLA-B*08:09', 'HLA-B*95:19', 'HLA-B*40:29', 'HLA-A*03:22', 'HLA-B*58:06', 'HLA-B*51:35',
                           'HLA-A*30:14L', 'HLA-B*39:16', 'HLA-B*15:89', 'HLA-B*57:02', 'HLA-A*34:05', 'HLA-B*15:06',
                           'HLA-B*44:30', 'HLA-A*11:29', 'HLA-A*11:26', 'HLA-B*40:10', 'HLA-B*18:12', 'HLA-B*15:68',
                           'HLA-B*40:47', 'HLA-A*68:07', 'HLA-B*54:02', 'HLA-B*44:32', 'HLA-B*42:04', 'HLA-B*39:05',
                           'HLA-B*35:31', 'HLA-B*40:16', 'HLA-B*18:14', 'HLA-B*44:26', 'HLA-C*03:18', 'HLA-C*04:01',
                           'HLA-B*07:08', 'HLA-B*08:10', 'HLA-B*07:47', 'HLA-B*52:03', 'HLA-C*03:31', 'HLA-C*01:04',
                           'HLA-A*02:62', 'HLA-B*51:06', 'HLA-B*95:15', 'HLA-B*35:72', 'HLA-B*51:22', 'HLA-C*04:08',
                           'HLA-C*06:04', 'HLA-A*68:21', 'HLA-C*01:06', 'HLA-A*11:10', 'HLA-A*24:46', 'HLA-B*95:17',
                           'HLA-A*26:30', 'HLA-C*02:09', 'HLA-B*55:16', 'HLA-B*39:18', 'HLA-C*05:13', 'HLA-A*23:05',
                           'HLA-A*11:12', 'HLA-B*47:03', 'HLA-C*05:09', 'HLA-A*26:32', 'HLA-A*02:48', 'HLA-B*44:38',
                           'HLA-B*51:45', 'HLA-B*07:45', 'HLA-A*02:68', 'HLA-B*37:11', 'HLA-C*03:12', 'HLA-A*30:19',
                           'HLA-A*32:14', 'HLA-B*15:28', 'HLA-A*11:20', 'HLA-A*11:06', 'HLA-A*68:23', 'HLA-A*03:26',
                           'HLA-B*58:02', 'HLA-B*07:31', 'HLA-C*03:35', 'HLA-B*40:23', 'HLA-C*03:10', 'HLA-B*15:23',
                           'HLA-B*47:05', 'HLA-B*42:06', 'HLA-B*40:14', 'HLA-A*26:16', 'HLA-A*01:01', 'HLA-B*07:22',
                           'HLA-A*30:11', 'HLA-B*27:30', 'HLA-C*07:41', 'HLA-B*15:25', 'HLA-A*11:15', 'HLA-C*05:02',
                           'HLA-B*38:12', 'HLA-A*01:03', 'HLA-A*23:04', 'HLA-B*07:24', 'HLA-A*30:17', 'HLA-B*58:09',
                           'HLA-C*08:13', 'HLA-A*31:15', 'HLA-A*02:44', 'HLA-C*06:06', 'HLA-B*42:02', 'HLA-B*40:18',
                           'HLA-C*15:05', 'HLA-B*44:16', 'HLA-B*07:43', 'HLA-B*44:27', 'HLA-B*38:10', 'HLA-C*01:08',
                           'HLA-B*44:03', 'HLA-A*30:15', 'HLA-B*40:26', 'HLA-B*44:36', 'HLA-B*40:64', 'HLA-A*68:36',
                           'HLA-B*55:10', 'HLA-B*44:14', 'HLA-B*18:05', 'HLA-A*32:01', 'HLA-B*15:93', 'HLA-B*44:05',
                           'HLA-B*95:13', 'HLA-B*35:70', 'HLA-B*08:31', 'HLA-C*04:06', 'HLA-A*24:34', 'HLA-B*56:09',
                           'HLA-B*39:12', 'HLA-C*01:03', 'HLA-B*15:77', 'HLA-B*27:31', 'HLA-B*35:54', 'HLA-A*24:51',
                           'HLA-A*24:10', 'HLA-B*42:09', 'HLA-B*27:20', 'HLA-B*15:14', 'HLA-B*55:04', 'HLA-B*07:14',
                           'HLA-B*14:04', 'HLA-B*51:05', 'HLA-B*35:52', 'HLA-B*51:42', 'HLA-C*02:07', 'HLA-B*51:17',
                           'HLA-A*03:07', 'HLA-B*18:10', 'HLA-A*11:16', 'HLA-B*55:02', 'HLA-A*02:40', 'HLA-B*40:62',
                           'HLA-A*11:01', 'HLA-A*24:14', 'HLA-A*03:05', 'HLA-B*38:04', 'HLA-A*11:18', 'HLA-B*15:62',
                           'HLA-A*02:46', 'HLA-B*40:06', 'HLA-B*51:09', 'HLA-B*07:41', 'HLA-B*15:21', 'HLA-A*02:89',
                           'HLA-B*39:19', 'HLA-A*01:07', 'HLA-B*57:10', 'HLA-B*07:20', 'HLA-B*07:11', 'HLA-C*03:33',
                           'HLA-C*02:11', 'HLA-B*49:01', 'HLA-B*15:45', 'HLA-B*44:33', 'HLA-A*23:01', 'HLA-A*24:33',
                           'HLA-B*44:07', 'HLA-B*48:09', 'HLA-B*07:04', 'HLA-B*40:68', 'HLA-B*07:29', 'HLA-B*44:34',
                           'HLA-C*08:14', 'HLA-B*15:43', 'HLA-A*23:03', 'HLA-B*44:09', 'HLA-B*39:31', 'HLA-A*29:05',
                           'HLA-B*15:07', 'HLA-A*02:33', 'HLA-B*27:14', 'HLA-B*55:24', 'HLA-A*03:09', 'HLA-B*35:13',
                           'HLA-A*25:03', 'HLA-B*40:55', 'HLA-B*08:02', 'HLA-B*48:05', 'HLA-B*15:64', 'HLA-A*29:06',
                           'HLA-B*41:03', 'HLA-A*02:86', 'HLA-C*01:07', 'HLA-C*03:08', 'HLA-A*26:17', 'HLA-C*02:13',
                           'HLA-B*49:03', 'HLA-B*08:04', 'HLA-A*29:09', 'HLA-B*15:66', 'HLA-B*15:27', 'HLA-B*44:13',
                           'HLA-B*27:10', 'HLA-A*66:04', 'HLA-C*01:05', 'HLA-B*51:01', 'HLA-A*32:08', 'HLA-B*35:56',
                           'HLA-B*40:40', 'HLA-A*02:81', 'HLA-C*15:11', 'HLA-A*31:01', 'HLA-A*24:30', 'HLA-C*07:17',
                           'HLA-C*12:07', 'HLA-B*15:60', 'HLA-B*51:10', 'HLA-C*07:08', 'HLA-A*23:09', 'HLA-B*38:09',
                           'HLA-C*04:04', 'HLA-A*24:68', 'HLA-B*07:30', 'HLA-B*40:27', 'HLA-A*68:35', 'HLA-B*54:01',
                           'HLA-B*35:38', 'HLA-B*56:01', 'HLA-B*56:12', 'HLA-A*26:10', 'HLA-A*26:15', 'HLA-B*55:17',
                           'HLA-A*02:70', 'HLA-B*51:21', 'HLA-B*49:05', 'HLA-B*15:49', 'HLA-B*56:03', 'HLA-A*25:02',
                           'HLA-C*05:03', 'HLA-A*29:03', 'HLA-A*02:31', 'HLA-A*26:19', 'HLA-A*31:03', 'HLA-A*26:13',
                           'HLA-C*02:17', 'HLA-B*15:47', 'HLA-B*39:35', 'HLA-B*51:15', 'HLA-B*35:59', 'HLA-C*05:05',
                           'HLA-A*68:22', 'HLA-A*29:01', 'HLA-B*39:11', 'HLA-B*07:06', 'HLA-A*68:32', 'HLA-B*73:01',
                           'HLA-B*07:27', 'HLA-B*48:01', 'HLA-B*52:05', 'HLA-A*01:19', 'HLA-C*04:23', 'HLA-A*66:06',
                           'HLA-A*23:06', 'HLA-B*35:16', 'HLA-A*31:08', 'HLA-B*40:42', 'HLA-B*44:47', 'HLA-B*55:23',
                           'HLA-C*05:08', 'HLA-C*04:21', 'HLA-B*39:17', 'HLA-B*35:14', 'HLA-B*39:33', 'HLA-B*08:13',
                           'HLA-C*12:02', 'HLA-B*15:08', 'HLA-B*51:12', 'HLA-A*03:23', 'HLA-B*08:06', 'HLA-B*18:15',
                           'HLA-B*39:15', 'HLA-C*15:15', 'HLA-A*32:09', 'HLA-B*40:46', 'HLA-A*24:49', 'HLA-C*12:03',
                           'HLA-B*40:20', 'HLA-A*11:28', 'HLA-B*39:39', 'HLA-B*48:03', 'HLA-B*45:01', 'HLA-C*06:07',
                           'HLA-B*39:13', 'HLA-C*15:17', 'HLA-B*15:88', 'HLA-A*36:01', 'HLA-A*24:32', 'HLA-B*07:09',
                           'HLA-B*44:49', 'HLA-B*95:03', 'HLA-A*02:50', 'HLA-B*27:29', 'HLA-B*44:31', 'HLA-C*07:23',
                           'HLA-B*35:51', 'HLA-A*26:31', 'HLA-C*07:28', 'HLA-A*11:24', 'HLA-A*02:85', 'HLA-B*27:33',
                           'HLA-B*56:16', 'HLA-B*07:46', 'HLA-B*40:69', 'HLA-B*08:15', 'HLA-A*68:30', 'HLA-C*07:03',
                           'HLA-B*35:35', 'HLA-B*46:03', 'HLA-B*07:05', 'HLA-A*03:25', 'HLA-C*06:14', 'HLA-A*11:03',
                           'HLA-B*56:07', 'HLA-C*03:16', 'HLA-C*07:01', 'HLA-C*15:09', 'HLA-B*40:67', 'HLA-B*18:04',
                           'HLA-B*44:43', 'HLA-B*55:21', 'HLA-C*07:42', 'HLA-B*57:11', 'HLA-B*40:45', 'HLA-B*44:39',
                           'HLA-C*18:01', 'HLA-A*68:34', 'HLA-C*07:07', 'HLA-B*48:04', 'HLA-A*11:13', 'HLA-A*23:02',
                           'HLA-B*35:12', 'HLA-A*26:33', 'HLA-A*29:07', 'HLA-B*54:04', 'HLA-B*07:48', 'HLA-B*35:06',
                           'HLA-C*07:05', 'HLA-B*35:33', 'HLA-B*37:10', 'HLA-B*40:35', 'HLA-B*35:10', 'HLA-B*15:29',
                           'HLA-B*46:09', 'HLA-C*15:13', 'HLA-A*24:18', 'HLA-A*02:24', 'HLA-A*68:05', 'HLA-B*44:10',
                           'HLA-C*04:24', 'HLA-B*18:19', 'HLA-A*26:09', 'HLA-B*40:24', 'HLA-B*59:02', 'HLA-B*40:60',
                           'HLA-A*25:06', 'HLA-A*24:53', 'HLA-A*24:42', 'HLA-A*43:01', 'HLA-C*08:10', 'HLA-B*44:45',
                           'HLA-C*06:12', 'HLA-A*02:08', 'HLA-B*07:18', 'HLA-B*51:04', 'HLA-A*26:34', 'HLA-B*83:01',
                           'HLA-A*02:11', 'HLA-A*01:10', 'HLA-A*32:05', 'HLA-B*27:23', 'HLA-C*01:09', 'HLA-B*15:20',
                           'HLA-B*44:02', 'HLA-A*02:06', 'HLA-B*15:90', 'HLA-C*04:07', 'HLA-B*15:24', 'HLA-A*02:13',
                           'HLA-B*51:40', 'HLA-B*07:42', 'HLA-A*02:26', 'HLA-A*01:02', 'HLA-B*44:41', 'HLA-B*44:04',
                           'HLA-A*02:04', 'HLA-A*24:19', 'HLA-B*35:71', 'HLA-B*56:08', 'HLA-A*36:04', 'HLA-C*17:01',
                           'HLA-A*24:35', 'HLA-B*07:40', 'HLA-B*40:19', 'HLA-B*44:17'])
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
    def supportedAlleles(self):
        """
        A list of supported :class:`~Fred2.Core.Allele.Allele` models
        """
        return self.__alleles

    @property
    def supportedLength(self):
        """
        A list of supported :class:`~Fred2.Core.Peptide.Peptide` lengths
        """
        return self.__supported_length

    def convert_alleles(self, alleles):
        """
        Converts :class:`~Fred2.Core.Allele.Allele` into the internal :class:`~Fred2.Core.Allele.Allele`
        representation of the predictor and returns a string representation

        :param alleles: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is needed
        :type alleles: list(:class:`~Fred2.Core.Allele.Allele`)
        :return: Returns a string representation of the input :class:`~Fred2.Core.Allele.Allele`
        :rtype: list(str)
        """
        return ["%s_%s%s" % (a.locus, a.supertype, a.subtype) for a in alleles]

    def encode(self, peptides, allele):
        """
        Encodes the input with binary sparse encoding of the :class:`~Fred2.Core.Peptide.Peptide`

        :param str peptides: A list of :class:`~Fred2.Core.Peptide.Peptide` sequences
        :param str allele: The HLA :class:`~Fred2.Core.Allele.Allele` represented by a string
        :return: Dictionary with :class:`~Fred2.Core.Peptide.Peptide` as key and feature encoding as value (see svmlight
                 encoding scheme http://svmlight.joachims.org/)
        :rtype: dict(:class:`~Fred2.Core.Peptide.Peptide`, (tuple(int, list(tuple(int,float))))`
        """
        pca = [{'A': 0.008, 'C': -0.132, 'E': 0.221, 'D': 0.303, 'G': 0.218, 'F': -0.329, 'I': -0.353, 'H': 0.023,
                'K': 0.243, 'M': -0.239, 'L': -0.267, 'N': 0.255, 'Q': 0.149, 'P': 0.173, 'S': 0.199, 'R': 0.171,
                'T': 0.068, 'W': -0.296, 'V': -0.274, 'Y': -0.141},
               {'A': 0.134, 'C': 0.174, 'E': -0.28, 'D': -0.057, 'G': 0.562, 'F': -0.023, 'I': 0.071, 'H': -0.177,
                'K': -0.339, 'M': -0.141, 'L': 0.018, 'N': 0.038, 'Q': -0.184, 'P': 0.286, 'S': 0.238, 'R': -0.361,
                'T': 0.147, 'W': -0.186, 'V': 0.136, 'Y': -0.057},
               {'A': -0.475, 'C': 0.07, 'E': -0.315, 'D': -0.014, 'G': -0.024, 'F': 0.072, 'I': -0.088, 'H': 0.041,
                'K': -0.044, 'M': -0.155, 'L': -0.265, 'N': 0.117, 'Q': -0.03, 'P': 0.407, 'S': -0.015, 'R': 0.107,
                'T': -0.015, 'W': 0.389, 'V': -0.187, 'Y': 0.425},
               {'A': -0.039, 'C': 0.565, 'E': 0.157, 'D': 0.225, 'G': 0.018, 'F': -0.002, 'I': -0.195, 'H': 0.28,
                'K': -0.325, 'M': 0.321, 'L': -0.274, 'N': 0.118, 'Q': 0.035, 'P': -0.215, 'S': -0.068, 'R': -0.258,
                'T': -0.132, 'W': 0.083, 'V': -0.196, 'Y': -0.096},
               {'A': 0.181, 'C': -0.374, 'E': 0.303, 'D': 0.156, 'G': 0.106, 'F': 0.208, 'I': -0.107, 'H': -0.021,
                'K': -0.027, 'M': 0.077, 'L': 0.206, 'N': -0.055, 'Q': -0.112, 'P': 0.384, 'S': -0.196, 'R': -0.364,
                'T': -0.274, 'W': 0.297, 'V': -0.299, 'Y': -0.091}]

        def __encode(pep, a):
            encoding = list(zip(list(range(1, 46)), UniTope_encodedAlleles[a + "_9"]))
            c = 46
            for p in str(pep):
                for i, pc in enumerate(pca):
                    encoding.append((c + i, pc[p]))
                c += 5
            return 0, encoding

        if isinstance(peptides, collections.Iterable):
            return {p: __encode(p, allele) for p in peptides}
        else:
            return {peptides: __encode(peptides, allele)}

    def predict(self, peptides, alleles=None, **kwargs):
        """
        Returns predictions for given peptides an alleles. If no alleles are given, predictions for all available models
        are made.

        :param peptides: A single :class:`~Fred2.Core.Peptide.Peptide` or a list of :class:`~Fred2.Core.Peptide.Peptide`
        :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`) or :class:`~Fred2.Core.Peptide.Peptide`
        :param alleles: A list of :class:`~Fred2.Core.Allele.Allele`
        :type alleles: list(:class:`~Fred2.Core.Allele.Allele`) or :class:`~Fred2.Core.Allele.Allele`
        :param kwargs: optional parameter (not used yet)
        :return: Returns a :class:`~Fred2.Core.Result.EpitopePredictionResult` object with the prediction results
        :rtype: :class:`~Fred2.Core.Result.EpitopePredictionResult`
        """
        if isinstance(peptides, Peptide):
            pep_seqs = {str(peptides): peptides}
        else:
            pep_seqs = {}
            for p in peptides:
                if not isinstance(p, Peptide):
                    raise ValueError("Input is not of type Protein or Peptide")
                pep_seqs[str(p)] = p

        chunksize = len(pep_seqs)
        if 'chunks' in kwargs:
            chunksize = kwargs['chunks']

        if alleles is None:
            al = [Allele(a) for a in self.supportedAlleles]
            allales_string = {conv_a: a for conv_a, a in zip(self.convert_alleles(al), al)}
        else:
            if isinstance(alleles, Allele):
                alleles = [alleles]
            if any(not isinstance(p, Allele) for p in alleles):
                raise ValueError("Input is not of type Allele")
            allales_string = {conv_a: a for conv_a, a in zip(self.convert_alleles(alleles), alleles)}

        # group peptides by length and
        result = {}

        model_path = pkg_resources.resource_filename("Fred2.Data.svms.%s" % self.name, "%s" % self.name)
        # model_path = os.path.abspath("../Data/svms/%s/%s"%(self.name, self.name))
        model = svmlight.read_model(model_path)

        for length, peps in itertools.groupby(iter(pep_seqs.keys()), key=lambda x: len(x)):
            # load svm model
            peps = list(peps)
            if length != 9:
                warnings.warn("Peptide length of %i is not supported by UniTope" % length)
                continue

            for a in list(allales_string.keys()):
                if allales_string[a].name in self.supportedAlleles:
                    for i in range(0, len(peps), chunksize):
                        encoding = self.encode(peps[i:i+chunksize], a)
                        pred = svmlight.classify(model, list(encoding.values()))
                        result[allales_string[a]] = {}
                        for pep, score in zip(list(encoding.keys()), pred):
                            result[allales_string[a]][pep_seqs[pep]] = score

        if not result:
            raise ValueError("No predictions could be made for given input. Check your \
            epitope length and HLA allele combination.")
        df_result = EpitopePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i, self.name)) for i in df_result.index],
                                                        names=['Seq', 'Method'])
        return df_result

# TODO: should we integrate this method or not? This means we have to HLA-DRag around ~500MB of model data just for this
# class MHCIIMulti(AEpitopePrediction, AExternal):
#     """
#         Implements MHCIIMulti
#     """
#
#     __name = "mhcIImulti"
#     __supported_length = frozenset([15])
#
#     models_dir = pkg_resources.resource_filename("Fred2.Data.svms.MHCIIMulti", "")
#     pockets = pkg_resources.resource_filename("Fred2.Data.svms.MHCIIMulti", "pockets.txt")
#     __command = "MHCIILeveraging %s %s %s "+str(models_dir)+" "+str(pockets)
#     __alleles = frozenset(['HLA-DRB3*02:21', 'HLA-DRB3*02:20', 'HLA-DRB1*14:22', 'HLA-DRB1*11:63', 'HLA-DRB1*11:62', 'HLA-DRB1*11:61', 'HLA-DRB1*11:60',
#                  'HLA-DRB1*11:67', 'HLA-DRB1*11:66', 'HLA-DRB1*11:64', 'HLA-DRB1*08:04:01', 'HLA-DRB1*08:04:02', 'HLA-DRB1*08:04:03',
#                  'HLA-DRB1*08:04:04', 'HLA-DRB1*04:59', 'HLA-DRB1*14:21', 'HLA-DRB1*04:54', 'HLA-DRB1*04:55', 'HLA-DRB1*04:56', 'HLA-DRB1*04:57',
#                  'HLA-DRB1*04:50', 'HLA-DRB1*04:51', 'HLA-DRB1*04:52', 'HLA-DRB1*04:53', 'HLA-DRB1*14:24', 'HLA-DRB1*14:25', 'HLA-DRB1*07:08',
#                  'HLA-DRB1*07:09', 'HLA-DRB1*13:29', 'HLA-DRB1*13:28', 'HLA-DRB1*13:25', 'HLA-DRB1*07:03', 'HLA-DRB1*13:27', 'HLA-DRB1*07:01',
#                  'HLA-DRB1*07:06', 'HLA-DRB1*07:07', 'HLA-DRB1*07:04', 'HLA-DRB1*07:05', 'HLA-DRB1*15:09', 'HLA-DRB1*01:02:04', 'HLA-DRB1*01:02:03',
#                  'HLA-DRB1*01:02:02', 'HLA-DRB1*01:02:01', 'HLA-DRB1*15:08', 'HLA-DRB1*03:37', 'HLA-DRB1*03:36', 'HLA-DRB1*03:35', 'HLA-DRB1*03:34',
#                  'HLA-DRB1*03:33', 'HLA-DRB1*03:32', 'HLA-DRB1*03:31', 'HLA-DRB1*03:30', 'HLA-DRB5*01:04', 'HLA-DRB5*01:05', 'HLA-DRB1*14:26',
#                  'HLA-DRB5*01:07', 'HLA-DRB5*01:01', 'HLA-DRB1*03:39', 'HLA-DRB1*03:38', 'HLA-DRB1*11:02:01', 'HLA-DRB1*11:02:02',
#                  'HLA-DRB3*01:01:02:02', 'HLA-DRB1*01:19', 'HLA-DRB1*01:18', 'HLA-DRB1*01:17', 'HLA-DRB1*01:16', 'HLA-DRB1*01:15', 'HLA-DRB1*01:14',
#                  'HLA-DRB1*01:13', 'HLA-DRB1*01:12', 'HLA-DRB1*01:11', 'HLA-DRB1*01:10', 'HLA-DRB3*02:17', 'HLA-DRB1*14:27', 'HLA-DRB1*11:51',
#                  'HLA-DRB1*11:12:01', 'HLA-DRB1*11:12:02', 'HLA-DRB1*14:59', 'HLA-DRB1*14:58', 'HLA-DRB1*14:51', 'HLA-DRB1*14:50', 'HLA-DRB1*14:53',
#                  'HLA-DRB1*14:52', 'HLA-DRB1*14:55', 'HLA-DRB1*14:54', 'HLA-DRB1*14:57', 'HLA-DRB1*14:56', 'HLA-DRB3*02:09', 'HLA-DRB3*02:08',
#                  'HLA-DRB3*02:03', 'HLA-DRB3*02:01', 'HLA-DRB3*02:07', 'HLA-DRB3*02:06', 'HLA-DRB3*02:05', 'HLA-DRB3*02:04', 'HLA-DRB1*09:01',
#                  'HLA-DRB1*09:02', 'HLA-DRB1*09:03', 'HLA-DRB1*09:04', 'HLA-DRB1*09:05', 'HLA-DRB1*09:06', 'HLA-DRB1*15:01:01', 'HLA-DRB1*15:01:03',
#                  'HLA-DRB1*15:01:02', 'HLA-DRB1*15:01:05', 'HLA-DRB1*15:01:04', 'HLA-DRB1*15:01:06', 'HLA-DRB1*08:32', 'HLA-DRB1*08:33',
#                  'HLA-DRB1*08:30', 'HLA-DRB1*08:31', 'HLA-DRB1*08:34', 'HLA-DRB1*04:18', 'HLA-DRB1*11:26', 'HLA-DRB1*11:25', 'HLA-DRB1*11:24',
#                  'HLA-DRB1*11:23', 'HLA-DRB1*11:22', 'HLA-DRB1*11:21', 'HLA-DRB1*11:20', 'HLA-DRB1*04:10', 'HLA-DRB1*04:11', 'HLA-DRB1*04:12',
#                  'HLA-DRB1*04:13', 'HLA-DRB1*04:14', 'HLA-DRB1*04:15', 'HLA-DRB1*11:29', 'HLA-DRB1*11:28', 'HLA-DRB3*02:10', 'HLA-DRB1*04:58',
#                  'HLA-DRB1*13:05:02', 'HLA-DRB1*13:05:01', 'HLA-DRB1*15:01', 'HLA-DRB1*15:03', 'HLA-DRB1*15:05', 'HLA-DRB1*15:04', 'HLA-DRB1*15:07',
#                  'HLA-DRB1*15:06', 'HLA-DRB1*04:01:03', 'HLA-DRB1*04:01:02', 'HLA-DRB1*04:01:01', 'HLA-DRB1*11:52', 'HLA-DRB1*11:53',
#                  'HLA-DRB1*11:50', 'HLA-DRB1*08:03:02', 'HLA-DRB1*11:56', 'HLA-DRB1*11:57', 'HLA-DRB1*11:55', 'HLA-DRB1*11:58', 'HLA-DRB1*11:59',
#                  'HLA-DRB1*13:82', 'HLA-DRB1*04:07:01', 'HLA-DRB1*04:07:03', 'HLA-DRB1*04:07:02', 'HLA-DRB1*13:81', 'HLA-DRB1*14:01:01',
#                  'HLA-DRB1*16:04', 'HLA-DRB1*12:13', 'HLA-DRB1*12:12', 'HLA-DRB1*12:11', 'HLA-DRB1*12:10', 'HLA-DRB1*03:02:01', 'HLA-DRB1*12:16',
#                  'HLA-DRB1*12:15', 'HLA-DRB1*03:02:02', 'HLA-DRB1*11:65:02', 'HLA-DRB1*11:65:01', 'HLA-DRB1*16:07', 'HLA-DRB1*13:58',
#                  'HLA-DRB1*13:59', 'HLA-DRB1*13:02:02', 'HLA-DRB1*13:02:03', 'HLA-DRB1*13:02:01', 'HLA-DRB1*11:08:02', 'HLA-DRB1*13:54',
#                  'HLA-DRB1*11:08:01', 'HLA-DRB1*03:08', 'HLA-DRB1*03:09', 'HLA-DRB1*03:06', 'HLA-DRB1*03:07', 'HLA-DRB1*03:04', 'HLA-DRB1*03:03',
#                  'HLA-DRB1*03:01', 'HLA-DRB1*14:15', 'HLA-DRB1*14:14', 'HLA-DRB1*04:05:04', 'HLA-DRB1*13:57', 'HLA-DRB1*14:11', 'HLA-DRB1*14:10',
#                  'HLA-DRB1*14:13', 'HLA-DRB1*14:12', 'HLA-DRB1*14:44:02', 'HLA-DRB1*14:44:01', 'HLA-DRB1*14:19', 'HLA-DRB1*14:18', 'HLA-DRB1*13:51',
#                  'HLA-DRB1*13:52', 'HLA-DRB1*12:14', 'HLA-DRB1*13:14:01', 'HLA-DRB1*13:14:02', 'HLA-DRB1*13:24', 'HLA-DRB1*13:26', 'HLA-DRB1*13:21',
#                  'HLA-DRB1*14:07:02', 'HLA-DRB1*13:20', 'HLA-DRB3*03:02', 'HLA-DRB3*03:03', 'HLA-DRB1*13:23', 'HLA-DRB1*13:22', 'HLA-DRB1*11:27:01',
#                  'HLA-DRB1*11:27:02', 'HLA-DRB1*14:05:01', 'HLA-DRB1*14:05:02', 'HLA-DRB1*14:05:03', 'HLA-DRB1*11:18', 'HLA-DRB1*15:17N',
#                  'HLA-DRB1*11:16', 'HLA-DRB1*11:17', 'HLA-DRB1*11:15', 'HLA-DRB1*11:13', 'HLA-DRB1*11:10', 'HLA-DRB1*04:43', 'HLA-DRB1*04:42',
#                  'HLA-DRB1*04:41', 'HLA-DRB1*04:40', 'HLA-DRB1*04:47', 'HLA-DRB1*04:46', 'HLA-DRB1*04:45', 'HLA-DRB1*04:44', 'HLA-DRB1*04:49',
#                  'HLA-DRB1*04:48', 'HLA-DRB1*10:02', 'HLA-DRB1*13:10', 'HLA-DRB1*13:11', 'HLA-DRB1*13:12', 'HLA-DRB1*13:13', 'HLA-DRB1*13:15',
#                  'HLA-DRB1*13:16', 'HLA-DRB1*13:17', 'HLA-DRB1*07:11', 'HLA-DRB1*13:19', 'HLA-DRB1*07:13', 'HLA-DRB1*07:12', 'HLA-DRB1*07:15',
#                  'HLA-DRB1*07:14', 'HLA-DRB4*01:03:01:01', 'HLA-DRB1*12:01:02', 'HLA-DRB1*16:01:02', 'HLA-DRB1*12:01:01', 'HLA-DRB1*16:01:01',
#                  'HLA-DRB5*01:13', 'HLA-DRB5*01:12', 'HLA-DRB5*01:11', 'HLA-DRB3*01:01:03', 'HLA-DRB3*01:01:04', 'HLA-DRB1*12:03:02',
#                  'HLA-DRB4*01:03:03', 'HLA-DRB4*01:03:02', 'HLA-DRB4*01:03:04', 'HLA-DRB1*03:05:01', 'HLA-DRB1*03:05:02', 'HLA-DRB1*13:80',
#                  'HLA-DRB1*14:16', 'HLA-DRB1*14:68', 'HLA-DRB1*14:69', 'HLA-DRB1*14:64', 'HLA-DRB1*14:65', 'HLA-DRB1*14:67', 'HLA-DRB1*14:60',
#                  'HLA-DRB1*14:61', 'HLA-DRB1*14:62', 'HLA-DRB1*14:63', 'HLA-DRB1*13:55', 'HLA-DRB5*01:06', 'HLA-DRB1*11:14:02', 'HLA-DRB1*11:14:01',
#                  'HLA-DRB1*13:56', 'HLA-DRB1*08:21', 'HLA-DRB1*08:20', 'HLA-DRB1*08:23', 'HLA-DRB1*08:22', 'HLA-DRB1*08:25', 'HLA-DRB1*08:24',
#                  'HLA-DRB1*08:27', 'HLA-DRB1*08:26', 'HLA-DRB1*08:29', 'HLA-DRB1*08:28', 'HLA-DRB5*01:03', 'HLA-DRB1*09:01:02', 'HLA-DRB1*09:01:03',
#                  'HLA-DRB5*01:09', 'HLA-DRB1*04:06:02', 'HLA-DRB1*04:09', 'HLA-DRB1*04:08', 'HLA-DRB1*04:05', 'HLA-DRB1*04:04', 'HLA-DRB1*04:02',
#                  'HLA-DRB1*04:01', 'HLA-DRB1*14:17', 'HLA-DRB3*01:01:02:01', 'HLA-DRB1*13:70', 'HLA-DRB1*11:41', 'HLA-DRB1*11:40', 'HLA-DRB1*11:43',
#                  'HLA-DRB1*11:42', 'HLA-DRB1*11:45', 'HLA-DRB1*11:44', 'HLA-DRB1*11:47', 'HLA-DRB1*11:46', 'HLA-DRB1*11:49', 'HLA-DRB1*11:48',
#                  'HLA-DRB1*04:70', 'HLA-DRB1*04:71', 'HLA-DRB1*04:03:04', 'HLA-DRB1*04:03:01', 'HLA-DRB1*04:03:03', 'HLA-DRB1*04:03:02',
#                  'HLA-DRB5*01:01:02', 'HLA-DRB5*01:01:01', 'HLA-DRB1*11:11:01', 'HLA-DRB1*11:11:02', 'HLA-DRB1*13:48', 'HLA-DRB1*13:43',
#                  'HLA-DRB1*13:42', 'HLA-DRB1*13:41', 'HLA-DRB1*13:40', 'HLA-DRB1*13:47', 'HLA-DRB1*13:46', 'HLA-DRB1*13:45', 'HLA-DRB1*13:44',
#                  'HLA-DRB1*13:53', 'HLA-DRB1*14:03:02', 'HLA-DRB1*14:03:01', 'HLA-DRB1*14:28', 'HLA-DRB1*14:29', 'HLA-DRB3*02:23', 'HLA-DRB3*02:22',
#                  'HLA-DRB1*03:19', 'HLA-DRB1*03:18', 'HLA-DRB1*03:15', 'HLA-DRB1*03:14', 'HLA-DRB1*03:17', 'HLA-DRB1*03:16', 'HLA-DRB1*03:11',
#                  'HLA-DRB1*03:10', 'HLA-DRB1*03:13', 'HLA-DRB1*03:12', 'HLA-DRB3*02:14', 'HLA-DRB1*16:03', 'HLA-DRB1*08:02:02', 'HLA-DRB1*08:02:03',
#                  'HLA-DRB1*08:02:01', 'HLA-DRB1*16:08', 'HLA-DRB1*16:09', 'HLA-DRB3*02:15', 'HLA-DRB1*13:65', 'HLA-DRB1*13:64', 'HLA-DRB1*15:02:04',
#                  'HLA-DRB1*15:02:01', 'HLA-DRB1*15:02:02', 'HLA-DRB1*15:02:03', 'HLA-DRB1*04:06:01', 'HLA-DRB5*01:02', 'HLA-DRB1*13:69',
#                  'HLA-DRB3*02:16', 'HLA-DRB1*13:68', 'HLA-DRB1*13:01:03', 'HLA-DRB1*13:01:02', 'HLA-DRB1*13:01:01', 'HLA-DRB1*11:09',
#                  'HLA-DRB1*11:05', 'HLA-DRB1*11:07', 'HLA-DRB1*11:01', 'HLA-DRB1*11:03', 'HLA-DRB1*04:36', 'HLA-DRB1*04:37', 'HLA-DRB1*04:34',
#                  'HLA-DRB1*04:35', 'HLA-DRB1*04:32', 'HLA-DRB1*04:33', 'HLA-DRB1*04:30', 'HLA-DRB1*04:31', 'HLA-DRB1*04:38', 'HLA-DRB1*04:39',
#                  'HLA-DRB3*02:11', 'HLA-DRB1*13:06', 'HLA-DRB1*13:04', 'HLA-DRB1*13:02', 'HLA-DRB1*13:01', 'HLA-DRB1*13:09', 'HLA-DRB1*13:08',
#                  'HLA-DRB3*02:12', 'HLA-DRB1*16:02:02', 'HLA-DRB1*16:02:01', 'HLA-DRB3*01:06', 'HLA-DRB3*01:07', 'HLA-DRB3*01:04', 'HLA-DRB3*01:05',
#                  'HLA-DRB3*01:02', 'HLA-DRB3*01:03', 'HLA-DRB3*01:01', 'HLA-DRB3*02:13', 'HLA-DRB3*01:08', 'HLA-DRB3*01:09', 'HLA-DRB1*14:07:01',
#                  'HLA-DRB1*13:72', 'HLA-DRB1*13:73', 'HLA-DRB1*11:04:04', 'HLA-DRB1*13:71', 'HLA-DRB1*11:04:02', 'HLA-DRB1*11:04:03',
#                  'HLA-DRB1*13:74', 'HLA-DRB1*11:04:01', 'HLA-DRB1*13:78', 'HLA-DRB1*13:79', 'HLA-DRB1*14:39', 'HLA-DRB1*13:07:01',
#                  'HLA-DRB1*13:07:02', 'HLA-DRB1*07:01:01', 'HLA-DRB1*07:01:02', 'HLA-DRB1*04:72', 'HLA-DRB1*14:73', 'HLA-DRB1*14:72',
#                  'HLA-DRB1*14:71', 'HLA-DRB1*14:70', 'HLA-DRB1*14:75', 'HLA-DRB1*14:74', 'HLA-DRB1*12:02:02', 'HLA-DRB1*12:02:01',
#                  'HLA-DRB1*11:01:05', 'HLA-DRB1*11:01:04', 'HLA-DRB1*11:01:07', 'HLA-DRB1*11:01:06', 'HLA-DRB1*11:01:01', 'HLA-DRB1*11:01:03',
#                  'HLA-DRB1*11:01:02', 'HLA-DRB1*04:05:05', 'HLA-DRB1*08:14', 'HLA-DRB1*08:15', 'HLA-DRB1*08:16', 'HLA-DRB1*08:17', 'HLA-DRB1*08:10',
#                  'HLA-DRB1*08:11', 'HLA-DRB1*08:12', 'HLA-DRB1*08:13', 'HLA-DRB1*04:05:03', 'HLA-DRB1*14:01:02', 'HLA-DRB1*14:01:03',
#                  'HLA-DRB1*08:18', 'HLA-DRB1*08:19', 'HLA-DRB1*04:05:02', 'HLA-DRB1*04:05:01', 'HLA-DRB1*11:19:01', 'HLA-DRB1*11:19:02',
#                  'HLA-DRB1*14:38', 'HLA-DRB1*14:37', 'HLA-DRB1*14:36', 'HLA-DRB3*02:18', 'HLA-DRB3*03:01:02', 'HLA-DRB3*03:01:01', 'HLA-DRB1*15:23',
#                  'HLA-DRB1*15:22', 'HLA-DRB1*15:21', 'HLA-DRB1*15:20', 'HLA-DRB1*15:27', 'HLA-DRB1*15:26', 'HLA-DRB1*15:25', 'HLA-DRB1*15:24',
#                  'HLA-DRB4*01:01:01:01', 'HLA-DRB1*04:69', 'HLA-DRB1*04:68', 'HLA-DRB1*04:61', 'HLA-DRB1*04:60', 'HLA-DRB1*04:63', 'HLA-DRB1*04:62',
#                  'HLA-DRB1*04:65', 'HLA-DRB1*04:64', 'HLA-DRB1*04:67', 'HLA-DRB1*04:66', 'HLA-DRB1*13:38', 'HLA-DRB1*13:39', 'HLA-DRB1*13:36',
#                  'HLA-DRB1*13:37', 'HLA-DRB1*13:34', 'HLA-DRB1*13:35', 'HLA-DRB1*13:32', 'HLA-DRB1*13:33', 'HLA-DRB1*13:30', 'HLA-DRB1*13:31',
#                  'HLA-DRB1*08:01:03', 'HLA-DRB1*08:01:02', 'HLA-DRB1*08:01:01', 'HLA-DRB1*13:50:01', 'HLA-DRB4*01:05', 'HLA-DRB1*03:01:01',
#                  'HLA-DRB1*03:01:02', 'HLA-DRB1*03:01:03', 'HLA-DRB1*03:01:04', 'HLA-DRB1*03:01:05', 'HLA-DRB1*03:01:06', 'HLA-DRB1*13:50:02',
#                  'HLA-DRB1*03:20', 'HLA-DRB1*03:21', 'HLA-DRB1*03:22', 'HLA-DRB1*03:23', 'HLA-DRB1*03:24', 'HLA-DRB1*03:25', 'HLA-DRB1*03:26',
#                  'HLA-DRB1*03:27', 'HLA-DRB1*03:28', 'HLA-DRB1*03:29', 'HLA-DRB1*14:35', 'HLA-DRB1*14:34', 'HLA-DRB1*14:33', 'HLA-DRB3*02:19',
#                  'HLA-DRB1*14:31', 'HLA-DRB1*14:30', 'HLA-DRB1*01:08', 'HLA-DRB1*01:09', 'HLA-DRB4*01:06', 'HLA-DRB1*04:73', 'HLA-DRB1*01:01',
#                  'HLA-DRB1*01:03', 'HLA-DRB1*01:04', 'HLA-DRB1*01:05', 'HLA-DRB1*01:06', 'HLA-DRB1*01:07', 'HLA-DRB1*16:05:02', 'HLA-DRB4*01:04',
#                  'HLA-DRB4*01:07', 'HLA-DRB1*16:05:01', 'HLA-DRB4*01:01', 'HLA-DRB4*01:02', 'HLA-DRB1*16:12', 'HLA-DRB1*16:11', 'HLA-DRB1*16:10',
#                  'HLA-DRB1*14:48', 'HLA-DRB1*14:49', 'HLA-DRB1*01:01:02', 'HLA-DRB1*01:01:03', 'HLA-DRB1*01:01:01', 'HLA-DRB1*14:42',
#                  'HLA-DRB1*14:43', 'HLA-DRB1*14:40', 'HLA-DRB1*14:41', 'HLA-DRB1*14:46', 'HLA-DRB1*14:47', 'HLA-DRB1*14:45', 'HLA-DRB1*11:54:02',
#                  'HLA-DRB1*11:54:01', 'HLA-DRB1*04:19', 'HLA-DRB1*13:03:01', 'HLA-DRB1*11:06:01', 'HLA-DRB1*11:06:02', 'HLA-DRB1*13:03:02',
#                  'HLA-DRB1*13:76', 'HLA-DRB1*13:77', 'HLA-DRB1*11:30', 'HLA-DRB1*11:31', 'HLA-DRB1*11:32', 'HLA-DRB1*11:33', 'HLA-DRB1*11:34',
#                  'HLA-DRB1*11:35', 'HLA-DRB1*11:36', 'HLA-DRB1*11:37', 'HLA-DRB1*11:38', 'HLA-DRB1*11:39', 'HLA-DRB1*13:75', 'HLA-DRB1*13:49',
#                  'HLA-DRB1*04:25', 'HLA-DRB1*04:24', 'HLA-DRB1*04:27', 'HLA-DRB1*04:26', 'HLA-DRB1*04:21', 'HLA-DRB1*04:20', 'HLA-DRB1*04:23',
#                  'HLA-DRB1*04:22', 'HLA-DRB1*04:29', 'HLA-DRB1*04:28', 'HLA-DRB1*14:32:02', 'HLA-DRB1*14:32:01', 'HLA-DRB1*04:16', 'HLA-DRB1*04:17',
#                  'HLA-DRB1*15:12', 'HLA-DRB1*15:13', 'HLA-DRB1*15:10', 'HLA-DRB1*15:11', 'HLA-DRB1*15:16', 'HLA-DRB1*15:14', 'HLA-DRB1*15:15',
#                  'HLA-DRB1*15:18', 'HLA-DRB1*15:19', 'HLA-DRB3*02:02:04', 'HLA-DRB3*02:02:05', 'HLA-DRB3*02:02:02', 'HLA-DRB3*02:02:03',
#                  'HLA-DRB3*02:02:01', 'HLA-DRB3*01:11', 'HLA-DRB3*01:10', 'HLA-DRB5*02:05', 'HLA-DRB5*02:04', 'HLA-DRB5*02:03', 'HLA-DRB5*02:02',
#                  'HLA-DRB1*13:61', 'HLA-DRB1*13:60', 'HLA-DRB1*13:63', 'HLA-DRB1*13:62', 'HLA-DRB1*12:08', 'HLA-DRB1*12:09', 'HLA-DRB1*13:67',
#                  'HLA-DRB1*13:66', 'HLA-DRB1*12:04', 'HLA-DRB1*12:05', 'HLA-DRB1*12:06', 'HLA-DRB1*12:07', 'HLA-DRB1*12:01', 'HLA-DRB5*01:08N',
#                  'HLA-DRB1*14:06', 'HLA-DRB1*13:18', 'HLA-DRB1*14:04', 'HLA-DRB1*14:02', 'HLA-DRB1*14:08', 'HLA-DRB1*14:09', 'HLA-DRB1*14:23:01',
#                  'HLA-DRB1*14:23:02', 'HLA-DRB1*08:09', 'HLA-DRB1*08:08', 'HLA-DRB1*10:01:01', 'HLA-DRB1*10:01:02', 'HLA-DRB1*08:02',
#                  'HLA-DRB1*08:01', 'HLA-DRB1*08:07', 'HLA-DRB1*08:06', 'HLA-DRB1*08:05', 'HLA-DRB1*14:20'])
#
#     @property
#     def name(self):
#         return self.__name
#
#     @property
#     def supportedAlleles(self):
#         return self.__alleles
#
#     @property
#     def supportedLength(self):
#         return self.__supported_length
#
#     @property
#     def command(self):
#         return self.__command
#
#     def parse_external_result(self, _file):
#         res = []
#         for l in _file:
#             if l.strip() == "":
#                 continue
#             res.append(float(l.strip()))
#         return res
#
#     def convert_alleles(self, alleles):
#         return ["%s_%s%s"%(a.locus, a.supertype, a.subtype) for a in alleles]
#
#     def predict(self, peptides, alleles=None, **kwargs):
#
#         if isinstance(peptides, Peptide):
#             pep_seqs = {str(peptides):peptides}
#         else:
#             if any(not isinstance(p, Peptide) for p in peptides):
#                 raise ValueError("Input is not of type Protein or Peptide")
#             pep_seqs = {str(p):p for p in peptides}
#
#         if alleles is None:
#             al = [Allele("HLA-"+a) for a in self.supportedAlleles]
#             allales_string = {conv_a:a for conv_a, a in itertools.izip(self.convert_alleles(al), al)}
#         else:
#             if isinstance(alleles, Allele):
#                 alleles = [alleles]
#             if any(not isinstance(p, Allele) for p in alleles):
#                 raise ValueError("Input is not of type Allele")
#             allales_string ={conv_a:a for conv_a, a in itertools.izip(self.convert_alleles(alleles),alleles)}
#
#         tmp_file = NamedTemporaryFile(delete=False)
#         pep = [ p for p in pep_seqs.keys() if len(p) >= max(self.__supported_length)]
#         if not pep:
#             raise ValueError("No epitopes with length >= %i"%max(self.__supported_length))
#
#         tmp_file.write("\n".join(pep_seqs.keys()))
#         tmp_file.close()
#         tmp_out = NamedTemporaryFile(delete=False)
#
#         results = {}
#         for a in allales_string.iterkeys():
#
#             #cmd = self.command%(data_file, a, prediction_file, self._modelpath) #modelpath?
#
#             r = subprocess.call(self.command%(tmp_file.name, a, tmp_out.name), shell=True)
#
#             if r == 127:
#                 raise RuntimeError("%s is not installed or globally executable."%self.name)
#             elif r == -6:
#                 warnings.warn("No model exists for allele %s."%str(allales_string[a]))
#                 continue
#             elif r != 0:
#                 warnings.warn("An unknown error occurred for method %s."%self.name)
#                 continue
#
#             results[allales_string[a]] = {p:s for p, s in itertools.izip(pep_seqs.values(), self.parse_external_result(tmp_out))}
#         print results
#         if any( not results[k] for k in results.iterkeys()):
#             raise ValueError("No predictions could be made for " +self.name+" given input. Check your "
#                              "epitope length and HLA allele combination.")
#         df_result = EpitopePredictionResult.from_dict(results)
#         df_result.index = pandas.MultiIndex.from_tuples([tuple((i,self.name)) for i in df_result.index],
#                                                         names=['Seq','Method'])
#         return df_result
