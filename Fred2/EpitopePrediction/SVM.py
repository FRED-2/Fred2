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

        :param peptides: A single :class:`~Fred2.Core.Peptide.Peptide` or a list of :class:`~Fred2.Core.Peptide.Peptide`s
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
            if any(not isinstance(p, Peptide) for p in peptides):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(p): p for p in peptides}

        if alleles is None:
            al = [Allele("HLA-" + a) for a in self.supportedAlleles]
            allales_string = {conv_a: a for conv_a, a in itertools.izip(self.convert_alleles(al), al)}
        else:
            if isinstance(alleles, Allele):
                alleles = [alleles]
            if any(not isinstance(p, Allele) for p in alleles):
                raise ValueError("Input is not of type Allele")
            allales_string = {conv_a: a for conv_a, a in itertools.izip(self.convert_alleles(alleles), alleles)}

        # group peptides by length and
        result = {}
        for length, peps in itertools.groupby(pep_seqs.iterkeys(), key=lambda x: len(x)):
            # load svm model

            if length not in self.supportedLength:
                warnings.warn("Peptide length of %i is not supported by %s" % (length, self.name))
                continue

            encoding = self.encode(peps)

            for a in allales_string.keys():
                model_path = pkg_resources.resource_filename("Fred2.Data.svms.%s" % self.name, "%s_%i" % (a, length))
                if not os.path.exists(model_path):
                    warnings.warn("No model exists for peptides of length %i or allele %s." % (length,
                                                                                               allales_string[a].name))
                    continue
                model = svmlight.read_model(model_path)

                model = svmlight.read_model(model_path)
                pred = svmlight.classify(model, encoding.values())
                result[allales_string[a]] = {}
                for pep, score in itertools.izip(encoding.keys(), pred):
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
    __alleles = frozenset(['A*02:01', 'A*02:01', 'A*11:01', 'A*11:01', 'A*24:02', 'B*15:01',
                           'B*15:01', 'B*18:01', 'B*18:01', 'B*27:05', 'B*35:01', 'B*37:01',
                           'B*51:01', 'B*51:01', 'C*04:01'])
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

        :param alleles: The :class:`~Fred2.Core.Allele.Allele`s for which the internal predictor representation is needed
        :type alleles: list(:class:`~Fred2.Core.Allele.Allele`)
        :return: Returns a string representation of the input :class:`~Fred2.Core.Allele.Allele`s
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
    __alleles = frozenset(['B*51:43', 'A*11:02', 'C*07:27', 'B*35:55', 'B*49:02', 'A*11:17',
                           'B*51:02', 'B*15:61', 'A*33:04', 'C*02:12', 'B*55:07', 'A*24:13',
                           'C*01:02', 'A*24:15', 'A*24:59', 'A*11:19', 'B*15:63', 'A*32:06',
                           'B*55:05', 'B*47:04', 'B*67:01', 'C*16:09', 'B*53:02', 'B*56:04',
                           'B*37:12', 'B*55:03', 'B*39:30', 'C*04:20', 'B*08:21', 'A*74:02',
                           'B*53:04', 'B*56:06', 'A*01:06', 'C*03:29', 'B*55:01', 'A*36:03',
                           'A*02:87', 'A*24:17', 'A*02:02', 'A*24:50', 'C*03:22', 'C*12:09',
                           'A*02:17', 'B*55:08', 'B*44:15', 'B*15:46', 'C*16:02', 'C*03:24',
                           'C*07:44', 'A*02:19', 'A*66:01', 'A*01:08', 'B*15:44', 'A*02:77',
                           'C*18:02', 'C*07:29', 'B*44:06', 'C*07:21', 'A*66:03', 'A*74:04',
                           'A*02:80', 'B*44:11', 'B*15:42', 'A*02:75', 'B*15:65', 'B*44:08',
                           'A*02:97', 'A*24:52', 'A*66:05', 'A*74:06', 'C*07:25', 'B*27:02',
                           'B*15:40', 'A*02:73', 'B*40:05', 'B*15:67', 'B*58:05', 'A*33:06',
                           'B*08:03', 'C*03:26', 'B*15:69', 'A*30:08', 'C*07:09', 'C*06:11',
                           'B*18:21', 'B*41:07', 'B*15:91', 'C*15:10', 'C*07:34', 'C*12:06',
                           'B*35:39', 'C*06:13', 'B*08:27', 'B*38:01', 'A*31:02', 'C*15:12',
                           'B*27:08', 'A*74:08', 'B*13:15', 'A*02:79', 'B*53:06', 'A*25:05',
                           'A*03:18', 'B*95:08', 'B*15:83', 'A*24:54', 'B*53:11', 'C*08:11',
                           'B*15:05', 'B*53:08', 'B*56:02', 'A*24:66', 'B*38:05', 'B*15:81',
                           'A*24:56', 'A*33:09', 'B*51:38', 'B*15:01', 'B*35:58', 'B*15:48',
                           'B*35:17', 'B*54:03', 'C*03:04', 'A*02:71', 'B*55:22', 'B*51:14',
                           'B*95:06', 'A*68:29', 'B*48:06', 'B*95:21', 'B*18:01', 'C*03:02',
                           'A*02:36', 'B*55:09', 'B*51:16', 'B*95:04', 'A*24:58', 'A*03:14',
                           'A*02:39', 'B*35:25', 'C*12:04', 'B*40:25', 'B*51:18', 'B*35:50',
                           'C*07:22', 'B*15:10', 'A*03:16', 'A*26:14', 'A*31:04', 'C*15:14',
                           'B*35:15', 'B*44:48', 'B*27:28', 'B*39:36', 'A*03:10', 'A*31:06',
                           'C*15:16', 'A*02:21', 'A*33:01', 'B*15:87', 'B*40:53', 'B*35:36',
                           'B*18:07', 'A*02:42', 'B*44:40', 'B*15:12', 'A*30:06', 'B*27:05',
                           'B*35:30', 'C*07:06', 'C*05:11', 'B*78:02', 'B*44:42', 'B*07:37',
                           'B*27:07', 'B*40:21', 'B*38:07', 'A*01:09', 'B*18:03', 'C*07:24',
                           'A*02:29', 'A*30:02', 'B*27:01', 'C*07:45', 'A*26:04', 'B*95:02',
                           'B*15:85', 'B*40:51', 'A*74:01', 'B*35:11', 'C*03:14', 'C*03:06',
                           'B*27:03', 'C*07:43', 'B*55:20', 'A*26:06', 'B*38:14', 'B*08:22',
                           'B*35:32', 'C*07:04', 'A*26:08', 'B*08:20', 'C*16:04', 'A*03:12',
                           'B*51:32', 'B*35:19', 'B*07:02', 'B*18:08', 'A*24:25', 'C*16:06',
                           'A*02:55', 'A*02:10', 'B*51:30', 'A*02:25', 'A*23:14', 'A*68:26',
                           'B*13:10', 'A*29:11', 'B*27:18', 'A*02:57', 'A*02:12', 'A*26:22',
                           'B*40:57', 'A*02:27', 'A*23:12', 'A*02:99', 'B*44:44', 'A*02:09',
                           'A*02:51', 'A*02:14', 'A*26:20', 'B*08:28', 'B*35:34', 'C*07:02',
                           'A*23:10', 'B*27:12', 'A*92:09', 'B*44:46', 'B*07:33', 'A*02:07',
                           'C*17:02', 'B*58:12', 'C*02:02', 'B*40:39', 'A*33:05', 'A*74:03',
                           'C*05:14', 'B*48:08', 'B*13:01', 'B*27:35', 'B*15:34', 'B*45:03',
                           'B*44:18', 'B*15:19', 'A*33:03', 'B*39:22', 'B*13:03', 'B*08:26',
                           'B*15:32', 'B*45:05', 'B*07:38', 'B*53:01', 'A*02:22', 'A*92:03',
                           'B*39:20', 'C*01:10', 'B*08:24', 'C*16:08', 'B*53:03', 'B*51:36',
                           'B*56:15', 'A*74:09', 'B*39:26', 'A*29:13', 'A*24:21', 'A*02:69',
                           'B*53:05', 'A*26:28', 'C*03:28', 'B*15:11', 'A*26:26', 'B*35:65Q',
                           'A*24:05', 'B*48:07', 'B*13:02', 'A*11:25', 'A*02:05', 'C*03:21',
                           'A*02:90', 'A*24:03', 'A*02:37', 'B*50:01', 'A*02:76', 'A*02:03',
                           'C*03:23', 'B*08:29', 'A*02:16', 'A*92:07', 'A*02:74', 'C*02:06',
                           'B*07:26', 'A*02:01', 'C*03:25', 'C*05:10', 'B*45:07', 'A*74:05',
                           'A*02:18', 'A*26:24', 'B*08:17', 'A*92:01', 'A*02:72', 'C*02:04',
                           'A*33:07', 'C*03:27', 'B*47:02', 'A*74:07', 'B*47:01', 'B*15:36',
                           'B*15:17', 'B*13:09', 'A*01:17', 'C*04:16', 'B*53:07', 'B*15:04',
                           'A*30:09', 'B*40:49', 'B*40:37', 'B*18:20', 'B*07:10', 'B*14:01',
                           'B*27:09', 'C*04:14', 'B*53:09', 'B*40:04', 'A*02:78', 'B*15:58',
                           'B*40:31', 'B*18:22', 'B*07:12', 'B*39:24', 'B*13:06', 'A*02:92',
                           'B*39:10', 'A*03:01', 'A*24:26', 'B*40:33', 'B*18:24', 'A*32:02',
                           'B*95:09', 'B*53:10', 'B*15:99', 'A*30:07', 'B*13:11', 'A*24:28',
                           'A*24:07', 'B*15:39', 'C*16:01', 'B*07:28', 'B*39:28', 'B*78:05',
                           'A*33:08', 'B*52:07', 'A*11:23', 'A*26:01', 'B*35:60', 'B*95:20',
                           'B*07:51', 'B*40:08', 'C*03:07', 'B*35:66', 'B*40:44', 'B*95:22',
                           'A*03:15', 'B*15:74', 'B*39:06', 'B*40:28', 'C*03:05', 'B*37:07',
                           'B*51:34', 'B*95:05', 'B*13:04', 'A*26:02', 'A*03:17', 'A*01:13',
                           'C*03:03', 'C*08:08', 'B*37:05', 'A*26:03', 'A*02:58', 'B*15:78',
                           'B*44:28', 'C*17:04', 'B*58:14', 'B*55:12', 'C*04:18', 'B*51:19',
                           'A*24:08', 'B*15:31', 'A*30:01', 'A*29:15', 'B*27:06', 'B*15:13',
                           'C*08:05', 'C*04:11', 'B*40:52', 'A*03:19', 'B*44:20', 'B*18:06',
                           'B*15:55', 'B*07:36', 'C*08:07', 'B*38:08', 'B*40:50', 'A*02:52',
                           'B*78:01', 'C*14:05', 'B*15:57', 'C*05:12', 'C*06:08', 'B*37:09',
                           'B*40:02', 'B*56:05', 'B*95:01', 'B*08:23', 'A*26:05', 'B*78:03',
                           'B*18:02', 'B*15:51', 'B*15:72', 'A*30:03', 'B*95:07', 'B*45:02',
                           'A*26:07', 'B*35:62', 'B*15:53', 'B*15:70', 'B*39:02', 'B*48:13',
                           'C*03:09', 'B*35:57', 'B*51:33', 'A*02:93', 'A*31:05', 'B*54:10',
                           'B*37:01', 'B*18:09', 'A*24:22', 'A*02:96', 'C*01:11', 'A*02:95',
                           'A*80:01', 'B*15:84', 'A*24:24', 'C*02:15', 'B*15:38', 'C*12:12',
                           'A*03:13', 'A*23:13', 'B*56:17', 'A*24:43', 'B*58:13', 'C*17:03',
                           'B*39:41', 'C*08:09', 'C*04:15', 'B*40:56', 'A*02:54', 'B*51:31',
                           'A*31:09', 'A*92:08', 'B*15:92', 'A*24:41', 'B*27:04', 'B*07:32',
                           'C*04:13', 'B*55:19', 'B*40:54', 'A*02:56', 'A*26:23', 'A*24:62',
                           'B*45:04', 'B*07:39', 'B*15:35', 'C*04:03', 'C*08:02', 'B*40:38',
                           'B*52:08', 'B*48:11', 'C*05:15', 'B*15:33', 'C*12:18', 'B*40:58',
                           'A*92:02', 'B*39:23', 'A*11:04', 'C*01:13', 'B*40:34', 'B*51:37',
                           'B*48:15', 'B*40:12', 'B*41:01', 'A*68:28', 'A*11:22', 'B*08:25',
                           'C*07:18', 'B*40:36', 'A*26:29', 'C*07:36', 'B*35:48', 'B*39:27',
                           'A*02:91', 'A*24:20', 'A*68:24', 'B*35:43', 'A*31:13', 'C*15:03',
                           'A*24:04', 'A*26:21', 'B*51:03', 'A*24:29', 'C*12:15', 'B*35:45',
                           'A*11:11', 'A*30:04', 'A*24:55', 'A*31:11', 'A*24:02', 'A*26:27',
                           'A*24:57', 'B*50:04', 'B*51:28', 'A*68:06', 'B*15:16', 'B*51:39',
                           'B*15:03', 'A*68:13', 'A*02:28', 'A*92:04', 'B*15:96', 'B*15:82',
                           'A*68:04', 'A*32:13', 'B*15:18', 'B*15:37', 'A*68:15', 'A*92:06',
                           'B*39:14', 'B*35:47', 'C*07:39', 'B*40:30', 'B*52:04', 'B*07:13',
                           'A*34:01', 'B*82:02', 'C*04:19', 'B*55:13', 'A*32:03', 'B*58:11',
                           'A*34:03', 'B*50:02', 'B*51:13', 'A*01:14', 'B*13:08', 'B*15:15',
                           'C*04:17', 'A*01:20', 'A*24:38', 'C*12:17', 'B*39:29', 'B*78:04',
                           'C*07:31', 'B*40:07', 'B*39:38Q', 'C*15:07', 'B*40:32', 'B*35:02',
                           'C*06:02', 'B*44:21', 'B*35:41', 'B*46:05', 'B*40:01', 'B*46:01',
                           'B*44:51', 'A*24:06', 'B*07:17', 'A*24:27', 'A*68:02', 'B*07:50',
                           'B*40:66', 'B*38:02', 'B*55:15', 'A*68:17', 'C*07:20', 'B*13:16',
                           'B*35:61', 'B*35:23', 'C*07:11', 'B*18:18', 'B*38:03', 'A*68:19',
                           'B*37:06', 'B*13:14', 'B*35:67', 'A*68:33', 'B*14:02', 'C*07:40',
                           'B*57:05', 'B*51:46', 'A*29:02', 'B*55:11', 'B*54:09', 'B*37:04',
                           'B*35:04', 'B*15:30', 'B*07:19', 'B*57:07', 'A*01:12', 'B*15:75',
                           'B*37:02', 'A*02:38', 'A*32:07', 'A*32:12', 'A*24:64', 'B*57:01',
                           'B*44:29', 'B*07:34', 'A*24:47', 'B*35:21', 'C*07:13', 'B*15:98',
                           'B*53:12', 'B*15:56', 'C*14:04', 'A*68:08', 'B*40:03', 'A*34:06',
                           'A*29:14', 'A*02:30', 'C*08:06', 'B*44:24', 'C*04:10', 'B*15:50',
                           'C*14:02', 'B*52:02', 'B*39:01', 'A*02:65', 'A*24:44', 'B*37:08',
                           'B*35:08', 'C*08:12', 'C*08:04', 'B*39:03', 'A*02:67', 'B*48:12',
                           'B*15:73', 'B*35:29', 'A*02:34', 'A*24:23', 'B*13:12', 'A*68:31',
                           'B*44:25', 'B*40:09', 'A*02:61', 'B*48:10', 'B*15:71', 'B*35:27',
                           'C*07:15', 'A*74:10', 'B*39:37', 'B*35:63', 'A*02:84', 'B*15:52',
                           'B*39:34', 'B*58:01', 'C*12:13', 'B*57:03', 'C*03:15', 'C*07:35',
                           'A*24:65', 'C*06:10', 'B*58:07', 'C*12:11', 'B*40:13', 'C*03:17',
                           'B*35:69', 'A*02:20', 'B*45:06', 'A*24:67', 'B*95:18', 'A*30:18',
                           'B*48:02', 'B*55:18', 'C*04:12', 'A*30:12', 'C*03:19', 'C*14:08',
                           'A*24:61', 'A*66:02', 'A*34:02', 'A*29:10', 'B*40:48', 'B*56:18',
                           'B*57:09', 'A*30:10', 'A*68:27', 'B*15:54', 'C*14:06', 'A*24:63',
                           'A*34:04', 'A*29:12', 'B*52:01', 'B*39:08', 'B*48:16', 'A*30:16',
                           'A*11:05', 'B*52:06', 'A*31:12', 'B*07:35', 'B*51:29', 'B*56:11',
                           'C*03:32', 'C*08:03', 'C*12:05', 'C*06:09', 'A*68:20', 'B*51:24',
                           'C*03:30', 'C*12:19', 'B*95:16', 'A*02:63', 'A*30:13', 'B*54:07',
                           'A*02:49', 'B*52:09', 'B*40:15', 'A*11:07', 'B*07:07', 'B*51:07',
                           'B*46:04', 'B*35:18', 'B*48:14', 'B*27:26', 'A*68:25', 'A*31:07',
                           'C*03:34', 'C*01:12', 'A*29:16', 'C*03:13', 'B*35:49', 'C*07:37',
                           'B*27:21', 'B*35:64', 'C*15:02', 'B*46:06', 'B*58:08', 'B*07:23',
                           'C*07:26', 'B*44:37', 'C*12:14', 'B*41:08', 'A*02:47', 'B*56:14',
                           'B*42:07', 'C*06:03', 'A*69:01', 'B*07:25', 'B*27:25', 'B*08:11',
                           'B*95:10', 'A*68:10', 'B*40:65', 'A*32:11Q', 'B*15:95', 'B*42:05',
                           'B*35:42', 'A*25:01', 'A*11:08', 'C*02:03', 'C*08:01', 'B*95:12',
                           'C*02:05', 'B*41:04', 'C*02:08', 'B*56:10', 'B*40:11', 'B*35:44',
                           'B*38:11', 'A*11:27', 'B*51:26', 'B*14:03', 'B*95:14', 'A*68:14',
                           'B*41:06', 'B*51:23', 'C*07:38', 'C*06:05', 'C*16:07', 'A*26:12',
                           'B*27:16', 'B*18:11', 'B*27:24', 'B*42:08', 'B*35:01', 'B*07:16',
                           'A*02:45', 'B*54:06', 'A*24:39', 'B*38:06', 'A*03:24', 'B*18:13',
                           'A*03:06', 'B*35:03', 'C*07:19', 'B*82:01', 'A*11:09', 'A*29:04',
                           'B*81:01', 'B*44:50', 'A*03:04', 'C*15:06', 'B*46:02', 'B*58:04',
                           'A*11:14', 'B*38:13', 'C*12:10', 'B*08:14', 'B*40:63', 'B*27:32',
                           'A*03:02', 'C*15:04', 'B*07:21', 'C*12:16', 'B*40:61', 'A*02:41',
                           'C*07:30', 'A*32:04', 'B*13:17', 'B*27:15', 'B*52:10', 'A*34:08',
                           'B*38:15', 'C*02:10', 'B*40:59', 'B*46:08', 'C*07:14', 'B*35:26',
                           'C*05:06', 'B*27:13', 'B*07:03', 'A*68:12', 'A*36:02', 'B*57:04',
                           'B*08:18', 'A*03:08', 'A*68:03', 'C*07:16', 'B*35:24', 'B*08:01',
                           'A*68:16', 'B*07:15', 'B*57:06', 'B*08:16', 'B*27:34', 'C*15:08',
                           'A*68:01', 'C*07:10', 'B*35:22', 'B*39:32', 'B*15:97', 'B*35:05',
                           'B*81:02', 'B*41:02', 'B*44:22', 'C*03:11', 'B*40:43', 'B*27:36',
                           'C*05:01', 'C*07:12', 'B*35:20', 'B*08:05', 'B*35:07', 'B*41:05',
                           'B*44:12', 'B*27:11', 'B*35:09', 'B*14:06', 'C*12:08', 'B*07:44',
                           'B*14:05', 'B*27:27', 'B*57:08', 'A*24:31', 'B*51:08', 'B*15:80',
                           'B*59:01', 'A*34:07', 'B*55:14', 'B*44:35', 'B*15:86', 'A*24:37',
                           'C*04:05', 'B*67:02', 'A*31:10', 'A*02:64', 'B*56:13', 'C*05:04',
                           'C*14:03', 'B*35:37', 'B*13:13', 'B*27:19', 'A*26:18', 'B*39:09',
                           'A*02:66', 'B*49:04', 'C*02:14', 'A*02:35', 'A*68:09', 'B*15:76',
                           'A*74:11', 'B*27:17', 'B*15:02', 'B*39:04', 'A*02:60', 'C*02:16',
                           'B*51:20', 'B*35:28', 'B*08:12', 'B*15:09', 'B*42:01', 'B*08:07',
                           'A*03:20', 'B*35:46', 'A*25:04', 'A*32:10', 'A*02:59', 'B*35:68',
                           'B*08:09', 'B*95:19', 'B*40:29', 'A*03:22', 'B*58:06', 'B*51:35',
                           'A*30:14L', 'B*39:16', 'B*15:89', 'B*57:02', 'A*34:05', 'B*15:06',
                           'B*44:30', 'A*11:29', 'A*11:26', 'B*40:10', 'B*18:12', 'B*15:68',
                           'B*40:47', 'A*68:07', 'B*54:02', 'B*44:32', 'B*42:04', 'B*39:05',
                           'B*35:31', 'B*40:16', 'B*18:14', 'B*44:26', 'C*03:18', 'C*04:01',
                           'B*07:08', 'B*08:10', 'B*07:47', 'B*52:03', 'C*03:31', 'C*01:04',
                           'A*02:62', 'B*51:06', 'B*95:15', 'B*35:72', 'B*51:22', 'C*04:08',
                           'C*06:04', 'A*68:21', 'C*01:06', 'A*11:10', 'A*24:46', 'B*95:17',
                           'A*26:30', 'C*02:09', 'B*55:16', 'B*39:18', 'C*05:13', 'A*23:05',
                           'A*11:12', 'B*47:03', 'C*05:09', 'A*26:32', 'A*02:48', 'B*44:38',
                           'B*51:45', 'B*07:45', 'A*02:68', 'B*37:11', 'C*03:12', 'A*30:19',
                           'A*32:14', 'B*15:28', 'A*11:20', 'A*11:06', 'A*68:23', 'A*03:26',
                           'B*58:02', 'B*07:31', 'C*03:35', 'B*40:23', 'C*03:10', 'B*15:23',
                           'B*47:05', 'B*42:06', 'B*40:14', 'A*26:16', 'A*01:01', 'B*07:22',
                           'A*30:11', 'B*27:30', 'C*07:41', 'B*15:25', 'A*11:15', 'C*05:02',
                           'B*38:12', 'A*01:03', 'A*23:04', 'B*07:24', 'A*30:17', 'B*58:09',
                           'C*08:13', 'A*31:15', 'A*02:44', 'C*06:06', 'B*42:02', 'B*40:18',
                           'C*15:05', 'B*44:16', 'B*07:43', 'B*44:27', 'B*38:10', 'C*01:08',
                           'B*44:03', 'A*30:15', 'B*40:26', 'B*44:36', 'B*40:64', 'A*68:36',
                           'B*55:10', 'B*44:14', 'B*18:05', 'A*32:01', 'B*15:93', 'B*44:05',
                           'B*95:13', 'B*35:70', 'B*08:31', 'C*04:06', 'A*24:34', 'B*56:09',
                           'B*39:12', 'C*01:03', 'B*15:77', 'B*27:31', 'B*35:54', 'A*24:51',
                           'A*24:10', 'B*42:09', 'B*27:20', 'B*15:14', 'B*55:04', 'B*07:14',
                           'B*14:04', 'B*51:05', 'B*35:52', 'B*51:42', 'C*02:07', 'B*51:17',
                           'A*03:07', 'B*18:10', 'A*11:16', 'B*55:02', 'A*02:40', 'B*40:62',
                           'A*11:01', 'A*24:14', 'A*03:05', 'B*38:04', 'A*11:18', 'B*15:62',
                           'A*02:46', 'B*40:06', 'B*51:09', 'B*07:41', 'B*15:21', 'A*02:89',
                           'B*39:19', 'A*01:07', 'B*57:10', 'B*07:20', 'B*07:11', 'C*03:33',
                           'C*02:11', 'B*49:01', 'B*15:45', 'B*44:33', 'A*23:01', 'A*24:33',
                           'B*44:07', 'B*48:09', 'B*07:04', 'B*40:68', 'B*07:29', 'B*44:34',
                           'C*08:14', 'B*15:43', 'A*23:03', 'B*44:09', 'B*39:31', 'A*29:05',
                           'B*15:07', 'A*02:33', 'B*27:14', 'B*55:24', 'A*03:09', 'B*35:13',
                           'A*25:03', 'B*40:55', 'B*08:02', 'B*48:05', 'B*15:64', 'A*29:06',
                           'B*41:03', 'A*02:86', 'C*01:07', 'C*03:08', 'A*26:17', 'C*02:13',
                           'B*49:03', 'B*08:04', 'A*29:09', 'B*15:66', 'B*15:27', 'B*44:13',
                           'B*27:10', 'A*66:04', 'C*01:05', 'B*51:01', 'A*32:08', 'B*35:56',
                           'B*40:40', 'A*02:81', 'C*15:11', 'A*31:01', 'A*24:30', 'C*07:17',
                           'C*12:07', 'B*15:60', 'B*51:10', 'C*07:08', 'A*23:09', 'B*38:09',
                           'C*04:04', 'A*24:68', 'B*07:30', 'B*40:27', 'A*68:35', 'B*54:01',
                           'B*35:38', 'B*56:01', 'B*56:12', 'A*26:10', 'A*26:15', 'B*55:17',
                           'A*02:70', 'B*51:21', 'B*49:05', 'B*15:49', 'B*56:03', 'A*25:02',
                           'C*05:03', 'A*29:03', 'A*02:31', 'A*26:19', 'A*31:03', 'A*26:13',
                           'C*02:17', 'B*15:47', 'B*39:35', 'B*51:15', 'B*35:59', 'C*05:05',
                           'A*68:22', 'A*29:01', 'B*39:11', 'B*07:06', 'A*68:32', 'B*73:01',
                           'B*07:27', 'B*48:01', 'B*52:05', 'A*01:19', 'C*04:23', 'A*66:06',
                           'A*23:06', 'B*35:16', 'A*31:08', 'B*40:42', 'B*44:47', 'B*55:23',
                           'C*05:08', 'C*04:21', 'B*39:17', 'B*35:14', 'B*39:33', 'B*08:13',
                           'C*12:02', 'B*15:08', 'B*51:12', 'A*03:23', 'B*08:06', 'B*18:15',
                           'B*39:15', 'C*15:15', 'A*32:09', 'B*40:46', 'A*24:49', 'C*12:03',
                           'B*40:20', 'A*11:28', 'B*39:39', 'B*48:03', 'B*45:01', 'C*06:07',
                           'B*39:13', 'C*15:17', 'B*15:88', 'A*36:01', 'A*24:32', 'B*07:09',
                           'B*44:49', 'B*95:03', 'A*02:50', 'B*27:29', 'B*44:31', 'C*07:23',
                           'B*35:51', 'A*26:31', 'C*07:28', 'A*11:24', 'A*02:85', 'B*27:33',
                           'B*56:16', 'B*07:46', 'B*40:69', 'B*08:15', 'A*68:30', 'C*07:03',
                           'B*35:35', 'B*46:03', 'B*07:05', 'A*03:25', 'C*06:14', 'A*11:03',
                           'B*56:07', 'C*03:16', 'C*07:01', 'C*15:09', 'B*40:67', 'B*18:04',
                           'B*44:43', 'B*55:21', 'C*07:42', 'B*57:11', 'B*40:45', 'B*44:39',
                           'C*18:01', 'A*68:34', 'C*07:07', 'B*48:04', 'A*11:13', 'A*23:02',
                           'B*35:12', 'A*26:33', 'A*29:07', 'B*54:04', 'B*07:48', 'B*35:06',
                           'C*07:05', 'B*35:33', 'B*37:10', 'B*40:35', 'B*35:10', 'B*15:29',
                           'B*46:09', 'C*15:13', 'A*24:18', 'A*02:24', 'A*68:05', 'B*44:10',
                           'C*04:24', 'B*18:19', 'A*26:09', 'B*40:24', 'B*59:02', 'B*40:60',
                           'A*25:06', 'A*24:53', 'A*24:42', 'A*43:01', 'C*08:10', 'B*44:45',
                           'C*06:12', 'A*02:08', 'B*07:18', 'B*51:04', 'A*26:34', 'B*83:01',
                           'A*02:11', 'A*01:10', 'A*32:05', 'B*27:23', 'C*01:09', 'B*15:20',
                           'B*44:02', 'A*02:06', 'B*15:90', 'C*04:07', 'B*15:24', 'A*02:13',
                           'B*51:40', 'B*07:42', 'A*02:26', 'A*01:02', 'B*44:41', 'B*44:04',
                           'A*02:04', 'A*24:19', 'B*35:71', 'B*56:08', 'A*36:04', 'C*17:01',
                           'A*24:35', 'B*07:40', 'B*40:19', 'B*44:17'])
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
        Converts :class:`~Fred2.Core.Allele.Allele`s into the internal :class:`~Fred2.Core.Allele.Allele`
        representation of the predictor and returns a string representation

        :param alleles: The :class:`~Fred2.Core.Allele.Allele`s for which the internal predictor representation is needed
        :type alleles: list(:class:`~Fred2.Core.Allele.Allele`)
        :return: Returns a string representation of the input :class:`~Fred2.Core.Allele.Allele`s
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
        :rtype: dict(:class:`~Fred2.Core.Peptide.Peptide`, (tuple(int, list(tuple(int,float))))
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
            encoding = zip(xrange(1, 46), UniTope_encodedAlleles[a + "_9"])
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

        :param peptides: A single :class:`~Fred2.Core.Peptide.Peptide` or a list of :class:`~Fred2.Core.Peptide.Peptide`s
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
            if any(not isinstance(p, Peptide) for p in peptides):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(p): p for p in peptides}

        if alleles is None:
            al = [Allele("HLA-" + a) for a in self.supportedAlleles]
            allales_string = {conv_a: a for conv_a, a in itertools.izip(self.convert_alleles(al), al)}
        else:
            if isinstance(alleles, Allele):
                alleles = [alleles]
            if any(not isinstance(p, Allele) for p in alleles):
                raise ValueError("Input is not of type Allele")
            allales_string = {conv_a: a for conv_a, a in itertools.izip(self.convert_alleles(alleles), alleles)}

        # group peptides by length and
        result = {}

        model_path = pkg_resources.resource_filename("Fred2.Data.svms.%s" % self.name, "%s" % self.name)
        # model_path = os.path.abspath("../Data/svms/%s/%s"%(self.name, self.name))
        model = svmlight.read_model(model_path)

        for length, peps in itertools.groupby(pep_seqs.iterkeys(), key=lambda x: len(x)):
            # load svm model
            peps = list(peps)
            if length != 9:
                warnings.warn("Peptide length of %i is not supported by UniTope" % length)
                continue

            for a in allales_string.keys():
                if allales_string[a].name in self.supportedAlleles:
                    encoding = self.encode(peps, a)
                    pred = svmlight.classify(model, encoding.values())
                    result[allales_string[a]] = {}
                    for pep, score in itertools.izip(encoding.keys(), pred):
                        result[allales_string[a]][pep_seqs[pep]] = score

        if not result:
            raise ValueError("No predictions could be made for given input. Check your \
            epitope length and HLA allele combination.")
        df_result = EpitopePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i, self.name)) for i in df_result.index],
                                                        names=['Seq', 'Method'])
        return df_result

# TODO: should we integrate this method or not? This means we have to drag around ~500MB of model data just for this
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
#     __alleles = frozenset(['DRB3*02:21', 'DRB3*02:20', 'DRB1*14:22', 'DRB1*11:63', 'DRB1*11:62', 'DRB1*11:61', 'DRB1*11:60',
#                  'DRB1*11:67', 'DRB1*11:66', 'DRB1*11:64', 'DRB1*08:04:01', 'DRB1*08:04:02', 'DRB1*08:04:03',
#                  'DRB1*08:04:04', 'DRB1*04:59', 'DRB1*14:21', 'DRB1*04:54', 'DRB1*04:55', 'DRB1*04:56', 'DRB1*04:57',
#                  'DRB1*04:50', 'DRB1*04:51', 'DRB1*04:52', 'DRB1*04:53', 'DRB1*14:24', 'DRB1*14:25', 'DRB1*07:08',
#                  'DRB1*07:09', 'DRB1*13:29', 'DRB1*13:28', 'DRB1*13:25', 'DRB1*07:03', 'DRB1*13:27', 'DRB1*07:01',
#                  'DRB1*07:06', 'DRB1*07:07', 'DRB1*07:04', 'DRB1*07:05', 'DRB1*15:09', 'DRB1*01:02:04', 'DRB1*01:02:03',
#                  'DRB1*01:02:02', 'DRB1*01:02:01', 'DRB1*15:08', 'DRB1*03:37', 'DRB1*03:36', 'DRB1*03:35', 'DRB1*03:34',
#                  'DRB1*03:33', 'DRB1*03:32', 'DRB1*03:31', 'DRB1*03:30', 'DRB5*01:04', 'DRB5*01:05', 'DRB1*14:26',
#                  'DRB5*01:07', 'DRB5*01:01', 'DRB1*03:39', 'DRB1*03:38', 'DRB1*11:02:01', 'DRB1*11:02:02',
#                  'DRB3*01:01:02:02', 'DRB1*01:19', 'DRB1*01:18', 'DRB1*01:17', 'DRB1*01:16', 'DRB1*01:15', 'DRB1*01:14',
#                  'DRB1*01:13', 'DRB1*01:12', 'DRB1*01:11', 'DRB1*01:10', 'DRB3*02:17', 'DRB1*14:27', 'DRB1*11:51',
#                  'DRB1*11:12:01', 'DRB1*11:12:02', 'DRB1*14:59', 'DRB1*14:58', 'DRB1*14:51', 'DRB1*14:50', 'DRB1*14:53',
#                  'DRB1*14:52', 'DRB1*14:55', 'DRB1*14:54', 'DRB1*14:57', 'DRB1*14:56', 'DRB3*02:09', 'DRB3*02:08',
#                  'DRB3*02:03', 'DRB3*02:01', 'DRB3*02:07', 'DRB3*02:06', 'DRB3*02:05', 'DRB3*02:04', 'DRB1*09:01',
#                  'DRB1*09:02', 'DRB1*09:03', 'DRB1*09:04', 'DRB1*09:05', 'DRB1*09:06', 'DRB1*15:01:01', 'DRB1*15:01:03',
#                  'DRB1*15:01:02', 'DRB1*15:01:05', 'DRB1*15:01:04', 'DRB1*15:01:06', 'DRB1*08:32', 'DRB1*08:33',
#                  'DRB1*08:30', 'DRB1*08:31', 'DRB1*08:34', 'DRB1*04:18', 'DRB1*11:26', 'DRB1*11:25', 'DRB1*11:24',
#                  'DRB1*11:23', 'DRB1*11:22', 'DRB1*11:21', 'DRB1*11:20', 'DRB1*04:10', 'DRB1*04:11', 'DRB1*04:12',
#                  'DRB1*04:13', 'DRB1*04:14', 'DRB1*04:15', 'DRB1*11:29', 'DRB1*11:28', 'DRB3*02:10', 'DRB1*04:58',
#                  'DRB1*13:05:02', 'DRB1*13:05:01', 'DRB1*15:01', 'DRB1*15:03', 'DRB1*15:05', 'DRB1*15:04', 'DRB1*15:07',
#                  'DRB1*15:06', 'DRB1*04:01:03', 'DRB1*04:01:02', 'DRB1*04:01:01', 'DRB1*11:52', 'DRB1*11:53',
#                  'DRB1*11:50', 'DRB1*08:03:02', 'DRB1*11:56', 'DRB1*11:57', 'DRB1*11:55', 'DRB1*11:58', 'DRB1*11:59',
#                  'DRB1*13:82', 'DRB1*04:07:01', 'DRB1*04:07:03', 'DRB1*04:07:02', 'DRB1*13:81', 'DRB1*14:01:01',
#                  'DRB1*16:04', 'DRB1*12:13', 'DRB1*12:12', 'DRB1*12:11', 'DRB1*12:10', 'DRB1*03:02:01', 'DRB1*12:16',
#                  'DRB1*12:15', 'DRB1*03:02:02', 'DRB1*11:65:02', 'DRB1*11:65:01', 'DRB1*16:07', 'DRB1*13:58',
#                  'DRB1*13:59', 'DRB1*13:02:02', 'DRB1*13:02:03', 'DRB1*13:02:01', 'DRB1*11:08:02', 'DRB1*13:54',
#                  'DRB1*11:08:01', 'DRB1*03:08', 'DRB1*03:09', 'DRB1*03:06', 'DRB1*03:07', 'DRB1*03:04', 'DRB1*03:03',
#                  'DRB1*03:01', 'DRB1*14:15', 'DRB1*14:14', 'DRB1*04:05:04', 'DRB1*13:57', 'DRB1*14:11', 'DRB1*14:10',
#                  'DRB1*14:13', 'DRB1*14:12', 'DRB1*14:44:02', 'DRB1*14:44:01', 'DRB1*14:19', 'DRB1*14:18', 'DRB1*13:51',
#                  'DRB1*13:52', 'DRB1*12:14', 'DRB1*13:14:01', 'DRB1*13:14:02', 'DRB1*13:24', 'DRB1*13:26', 'DRB1*13:21',
#                  'DRB1*14:07:02', 'DRB1*13:20', 'DRB3*03:02', 'DRB3*03:03', 'DRB1*13:23', 'DRB1*13:22', 'DRB1*11:27:01',
#                  'DRB1*11:27:02', 'DRB1*14:05:01', 'DRB1*14:05:02', 'DRB1*14:05:03', 'DRB1*11:18', 'DRB1*15:17N',
#                  'DRB1*11:16', 'DRB1*11:17', 'DRB1*11:15', 'DRB1*11:13', 'DRB1*11:10', 'DRB1*04:43', 'DRB1*04:42',
#                  'DRB1*04:41', 'DRB1*04:40', 'DRB1*04:47', 'DRB1*04:46', 'DRB1*04:45', 'DRB1*04:44', 'DRB1*04:49',
#                  'DRB1*04:48', 'DRB1*10:02', 'DRB1*13:10', 'DRB1*13:11', 'DRB1*13:12', 'DRB1*13:13', 'DRB1*13:15',
#                  'DRB1*13:16', 'DRB1*13:17', 'DRB1*07:11', 'DRB1*13:19', 'DRB1*07:13', 'DRB1*07:12', 'DRB1*07:15',
#                  'DRB1*07:14', 'DRB4*01:03:01:01', 'DRB1*12:01:02', 'DRB1*16:01:02', 'DRB1*12:01:01', 'DRB1*16:01:01',
#                  'DRB5*01:13', 'DRB5*01:12', 'DRB5*01:11', 'DRB3*01:01:03', 'DRB3*01:01:04', 'DRB1*12:03:02',
#                  'DRB4*01:03:03', 'DRB4*01:03:02', 'DRB4*01:03:04', 'DRB1*03:05:01', 'DRB1*03:05:02', 'DRB1*13:80',
#                  'DRB1*14:16', 'DRB1*14:68', 'DRB1*14:69', 'DRB1*14:64', 'DRB1*14:65', 'DRB1*14:67', 'DRB1*14:60',
#                  'DRB1*14:61', 'DRB1*14:62', 'DRB1*14:63', 'DRB1*13:55', 'DRB5*01:06', 'DRB1*11:14:02', 'DRB1*11:14:01',
#                  'DRB1*13:56', 'DRB1*08:21', 'DRB1*08:20', 'DRB1*08:23', 'DRB1*08:22', 'DRB1*08:25', 'DRB1*08:24',
#                  'DRB1*08:27', 'DRB1*08:26', 'DRB1*08:29', 'DRB1*08:28', 'DRB5*01:03', 'DRB1*09:01:02', 'DRB1*09:01:03',
#                  'DRB5*01:09', 'DRB1*04:06:02', 'DRB1*04:09', 'DRB1*04:08', 'DRB1*04:05', 'DRB1*04:04', 'DRB1*04:02',
#                  'DRB1*04:01', 'DRB1*14:17', 'DRB3*01:01:02:01', 'DRB1*13:70', 'DRB1*11:41', 'DRB1*11:40', 'DRB1*11:43',
#                  'DRB1*11:42', 'DRB1*11:45', 'DRB1*11:44', 'DRB1*11:47', 'DRB1*11:46', 'DRB1*11:49', 'DRB1*11:48',
#                  'DRB1*04:70', 'DRB1*04:71', 'DRB1*04:03:04', 'DRB1*04:03:01', 'DRB1*04:03:03', 'DRB1*04:03:02',
#                  'DRB5*01:01:02', 'DRB5*01:01:01', 'DRB1*11:11:01', 'DRB1*11:11:02', 'DRB1*13:48', 'DRB1*13:43',
#                  'DRB1*13:42', 'DRB1*13:41', 'DRB1*13:40', 'DRB1*13:47', 'DRB1*13:46', 'DRB1*13:45', 'DRB1*13:44',
#                  'DRB1*13:53', 'DRB1*14:03:02', 'DRB1*14:03:01', 'DRB1*14:28', 'DRB1*14:29', 'DRB3*02:23', 'DRB3*02:22',
#                  'DRB1*03:19', 'DRB1*03:18', 'DRB1*03:15', 'DRB1*03:14', 'DRB1*03:17', 'DRB1*03:16', 'DRB1*03:11',
#                  'DRB1*03:10', 'DRB1*03:13', 'DRB1*03:12', 'DRB3*02:14', 'DRB1*16:03', 'DRB1*08:02:02', 'DRB1*08:02:03',
#                  'DRB1*08:02:01', 'DRB1*16:08', 'DRB1*16:09', 'DRB3*02:15', 'DRB1*13:65', 'DRB1*13:64', 'DRB1*15:02:04',
#                  'DRB1*15:02:01', 'DRB1*15:02:02', 'DRB1*15:02:03', 'DRB1*04:06:01', 'DRB5*01:02', 'DRB1*13:69',
#                  'DRB3*02:16', 'DRB1*13:68', 'DRB1*13:01:03', 'DRB1*13:01:02', 'DRB1*13:01:01', 'DRB1*11:09',
#                  'DRB1*11:05', 'DRB1*11:07', 'DRB1*11:01', 'DRB1*11:03', 'DRB1*04:36', 'DRB1*04:37', 'DRB1*04:34',
#                  'DRB1*04:35', 'DRB1*04:32', 'DRB1*04:33', 'DRB1*04:30', 'DRB1*04:31', 'DRB1*04:38', 'DRB1*04:39',
#                  'DRB3*02:11', 'DRB1*13:06', 'DRB1*13:04', 'DRB1*13:02', 'DRB1*13:01', 'DRB1*13:09', 'DRB1*13:08',
#                  'DRB3*02:12', 'DRB1*16:02:02', 'DRB1*16:02:01', 'DRB3*01:06', 'DRB3*01:07', 'DRB3*01:04', 'DRB3*01:05',
#                  'DRB3*01:02', 'DRB3*01:03', 'DRB3*01:01', 'DRB3*02:13', 'DRB3*01:08', 'DRB3*01:09', 'DRB1*14:07:01',
#                  'DRB1*13:72', 'DRB1*13:73', 'DRB1*11:04:04', 'DRB1*13:71', 'DRB1*11:04:02', 'DRB1*11:04:03',
#                  'DRB1*13:74', 'DRB1*11:04:01', 'DRB1*13:78', 'DRB1*13:79', 'DRB1*14:39', 'DRB1*13:07:01',
#                  'DRB1*13:07:02', 'DRB1*07:01:01', 'DRB1*07:01:02', 'DRB1*04:72', 'DRB1*14:73', 'DRB1*14:72',
#                  'DRB1*14:71', 'DRB1*14:70', 'DRB1*14:75', 'DRB1*14:74', 'DRB1*12:02:02', 'DRB1*12:02:01',
#                  'DRB1*11:01:05', 'DRB1*11:01:04', 'DRB1*11:01:07', 'DRB1*11:01:06', 'DRB1*11:01:01', 'DRB1*11:01:03',
#                  'DRB1*11:01:02', 'DRB1*04:05:05', 'DRB1*08:14', 'DRB1*08:15', 'DRB1*08:16', 'DRB1*08:17', 'DRB1*08:10',
#                  'DRB1*08:11', 'DRB1*08:12', 'DRB1*08:13', 'DRB1*04:05:03', 'DRB1*14:01:02', 'DRB1*14:01:03',
#                  'DRB1*08:18', 'DRB1*08:19', 'DRB1*04:05:02', 'DRB1*04:05:01', 'DRB1*11:19:01', 'DRB1*11:19:02',
#                  'DRB1*14:38', 'DRB1*14:37', 'DRB1*14:36', 'DRB3*02:18', 'DRB3*03:01:02', 'DRB3*03:01:01', 'DRB1*15:23',
#                  'DRB1*15:22', 'DRB1*15:21', 'DRB1*15:20', 'DRB1*15:27', 'DRB1*15:26', 'DRB1*15:25', 'DRB1*15:24',
#                  'DRB4*01:01:01:01', 'DRB1*04:69', 'DRB1*04:68', 'DRB1*04:61', 'DRB1*04:60', 'DRB1*04:63', 'DRB1*04:62',
#                  'DRB1*04:65', 'DRB1*04:64', 'DRB1*04:67', 'DRB1*04:66', 'DRB1*13:38', 'DRB1*13:39', 'DRB1*13:36',
#                  'DRB1*13:37', 'DRB1*13:34', 'DRB1*13:35', 'DRB1*13:32', 'DRB1*13:33', 'DRB1*13:30', 'DRB1*13:31',
#                  'DRB1*08:01:03', 'DRB1*08:01:02', 'DRB1*08:01:01', 'DRB1*13:50:01', 'DRB4*01:05', 'DRB1*03:01:01',
#                  'DRB1*03:01:02', 'DRB1*03:01:03', 'DRB1*03:01:04', 'DRB1*03:01:05', 'DRB1*03:01:06', 'DRB1*13:50:02',
#                  'DRB1*03:20', 'DRB1*03:21', 'DRB1*03:22', 'DRB1*03:23', 'DRB1*03:24', 'DRB1*03:25', 'DRB1*03:26',
#                  'DRB1*03:27', 'DRB1*03:28', 'DRB1*03:29', 'DRB1*14:35', 'DRB1*14:34', 'DRB1*14:33', 'DRB3*02:19',
#                  'DRB1*14:31', 'DRB1*14:30', 'DRB1*01:08', 'DRB1*01:09', 'DRB4*01:06', 'DRB1*04:73', 'DRB1*01:01',
#                  'DRB1*01:03', 'DRB1*01:04', 'DRB1*01:05', 'DRB1*01:06', 'DRB1*01:07', 'DRB1*16:05:02', 'DRB4*01:04',
#                  'DRB4*01:07', 'DRB1*16:05:01', 'DRB4*01:01', 'DRB4*01:02', 'DRB1*16:12', 'DRB1*16:11', 'DRB1*16:10',
#                  'DRB1*14:48', 'DRB1*14:49', 'DRB1*01:01:02', 'DRB1*01:01:03', 'DRB1*01:01:01', 'DRB1*14:42',
#                  'DRB1*14:43', 'DRB1*14:40', 'DRB1*14:41', 'DRB1*14:46', 'DRB1*14:47', 'DRB1*14:45', 'DRB1*11:54:02',
#                  'DRB1*11:54:01', 'DRB1*04:19', 'DRB1*13:03:01', 'DRB1*11:06:01', 'DRB1*11:06:02', 'DRB1*13:03:02',
#                  'DRB1*13:76', 'DRB1*13:77', 'DRB1*11:30', 'DRB1*11:31', 'DRB1*11:32', 'DRB1*11:33', 'DRB1*11:34',
#                  'DRB1*11:35', 'DRB1*11:36', 'DRB1*11:37', 'DRB1*11:38', 'DRB1*11:39', 'DRB1*13:75', 'DRB1*13:49',
#                  'DRB1*04:25', 'DRB1*04:24', 'DRB1*04:27', 'DRB1*04:26', 'DRB1*04:21', 'DRB1*04:20', 'DRB1*04:23',
#                  'DRB1*04:22', 'DRB1*04:29', 'DRB1*04:28', 'DRB1*14:32:02', 'DRB1*14:32:01', 'DRB1*04:16', 'DRB1*04:17',
#                  'DRB1*15:12', 'DRB1*15:13', 'DRB1*15:10', 'DRB1*15:11', 'DRB1*15:16', 'DRB1*15:14', 'DRB1*15:15',
#                  'DRB1*15:18', 'DRB1*15:19', 'DRB3*02:02:04', 'DRB3*02:02:05', 'DRB3*02:02:02', 'DRB3*02:02:03',
#                  'DRB3*02:02:01', 'DRB3*01:11', 'DRB3*01:10', 'DRB5*02:05', 'DRB5*02:04', 'DRB5*02:03', 'DRB5*02:02',
#                  'DRB1*13:61', 'DRB1*13:60', 'DRB1*13:63', 'DRB1*13:62', 'DRB1*12:08', 'DRB1*12:09', 'DRB1*13:67',
#                  'DRB1*13:66', 'DRB1*12:04', 'DRB1*12:05', 'DRB1*12:06', 'DRB1*12:07', 'DRB1*12:01', 'DRB5*01:08N',
#                  'DRB1*14:06', 'DRB1*13:18', 'DRB1*14:04', 'DRB1*14:02', 'DRB1*14:08', 'DRB1*14:09', 'DRB1*14:23:01',
#                  'DRB1*14:23:02', 'DRB1*08:09', 'DRB1*08:08', 'DRB1*10:01:01', 'DRB1*10:01:02', 'DRB1*08:02',
#                  'DRB1*08:01', 'DRB1*08:07', 'DRB1*08:06', 'DRB1*08:05', 'DRB1*14:20'])
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
