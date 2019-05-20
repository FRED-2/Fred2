# coding=utf-8
# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: EpitopePrediction.ANN
   :synopsis: This module contains all classes for ANN-based epitope prediction.
.. moduleauthor:: heumos

"""
import subprocess
import csv
import tempfile
from abc import abstractmethod

import pandas

from mhcflurry import Class1AffinityPredictor
from mhcnuggets.src.predict import predict as mhcnuggets_predict
from Fred2.Core import EpitopePredictionResult
from Fred2.Core.Base import AEpitopePrediction

import inspect


class BadSignatureException(Exception):
    pass


class SignatureCheckerMeta(type):
    def __new__(cls, name, baseClasses, d):
        # For each method in d, check to see if any base class already
        # defined a method with that name. If so, make sure the
        # signatures are the same.
        for methodName in d:
            f = d[methodName]
            for baseClass in baseClasses:
                try:
                    fBase = getattr(baseClass, methodName).__func__
                    if not inspect.getargspec(f) == inspect.getargspec(fBase):
                        raise BadSignatureException(str(methodName))
                except AttributeError:
                    # This method was not defined in this base class,
                    # So just go to the next base class.
                    continue

        return type(name, baseClasses, d)


class AANNEpitopePrediction(AEpitopePrediction):
    """
        Abstract base class for ANN predictions.
        Implements predict functionality
    """

    @abstractmethod
    def predict(self, peptides, alleles=None, binary=False, **kwargs):
        """
        All ANN based predictors have to implement their custom predict method.
        Furthermore, all of them have to use the metaclass SignatureCheckerMeta to check for any contract violations.
        They have to adhere to the following contract

        :param peptides: A single :class:`~Fred2.Core.Peptide.Peptide` or a list of :class:`~Fred2.Core.Peptide.Peptide`
        :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`) or :class:`~Fred2.Core.Peptide.Peptide`
        :param kwargs: optional parameter (not used yet)
        :return: Returns a :class:`pandas.DataFrame` object with the prediction results
        :rtype: :class:`pandas.DataFrame`
        """


try:
    class MHCNuggetsPredictor_1(AANNEpitopePrediction):
        """
        Implements MHCNuggets Class I

        .. note::
            Evaluation of machine learning methods to predict peptide binding to MHC Class I proteins
            Rohit Bhattacharya, Ashok Sivakumar, Collin Tokheim, Violeta Beleva Guthrie, Valsamo Anagnostou,
            Victor E. Velculescu, Rachel Karchin (2017) bioRxiv
        """
        __metaclass__ = SignatureCheckerMeta
        __alleles = frozenset(
            ["HLA-A01:01", "HLA-A02:01", "HLA-A02:02", "HLA-A02:03", "HLA-A02:04", "HLA-A02:05", "HLA-A02:06",
             "HLA-A02:07", "HLA-A02:08", "HLA-A02:09", "HLA-A02:10", "HLA-A02:11", "HLA-A02:12", "HLA-A02:14",
             "HLA-A02:16", "HLA-A02:17", "HLA-A02:19", "HLA-A02:50", "HLA-A03:01", "HLA-A03:02", "HLA-A03:19",
             "HLA-A11:01", "HLA-A11:02", "HLA-A23:01", "HLA-A24:01", "HLA-A24:02", "HLA-A24:03", "HLA-A25:01",
             "HLA-A26:01", "HLA-A26:02", "HLA-A26:03", "HLA-A29:01", "HLA-A29:02", "HLA-A30:01", "HLA-A30:02",
             "HLA-A30:03", "HLA-A30:04", "HLA-A31:01", "HLA-A32:01", "HLA-A32:07", "HLA-A32:15", "HLA-A33:01",
             "HLA-A33:03", "HLA-A66:01", "HLA-A68:01", "HLA-A68:02", "HLA-A68:23", "HLA-A69:01", "HLA-A74:01",
             "HLA-A80:01", "HLA-B07:01", "HLA-B07:02", "HLA-B08:01", "HLA-B08:02", "HLA-B08:03", "HLA-B12:01",
             "HLA-B13:02", "HLA-B14:01", "HLA-B14:02", "HLA-B15:01", "HLA-B15:02", "HLA-B15:03", "HLA-B15:08",
             "HLA-B15:09", "HLA-B15:10", "HLA-B15:13", "HLA-B15:16", "HLA-B15:17", "HLA-B15:42", "HLA-B18:01",
             "HLA-B27:01", "HLA-B27:02", "HLA-B27:03", "HLA-B27:04", "HLA-B27:05", "HLA-B27:06", "HLA-B27:09",
             "HLA-B27:10", "HLA-B27:20", "HLA-B35:01", "HLA-B35:02", "HLA-B35:03", "HLA-B35:08", "HLA-B37:01",
             "HLA-B38:01", "HLA-B39:01", "HLA-B39:06", "HLA-B39:09", "HLA-B39:10", "HLA-B40:01", "HLA-B40:02",
             "HLA-B40:13", "HLA-B41:03", "HLA-B41:04", "HLA-B42:01", "HLA-B42:02", "HLA-B44:01", "HLA-B44:02",
             "HLA-B44:03", "HLA-B44:05", "HLA-B45:01", "HLA-B45:06", "HLA-B46:01", "HLA-B48:01", "HLA-B51:01",
             "HLA-B51:02", "HLA-B52:01", "HLA-B53:01", "HLA-B54:01", "HLA-B55:01", "HLA-B55:02", "HLA-B56:01",
             "HLA-B57:01", "HLA-B57:02", "HLA-B57:03", "HLA-B58:01", "HLA-B58:02", "HLA-B60:01", "HLA-B61:01",
             "HLA-B62:01", "HLA-B73:01", "HLA-B81:01", "HLA-B83:01"])
        __supported_length = frozenset([5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22])
        __name = "mhcnuggets-class-1"
        __version = "2.0"

        # the interface defines three class properties
        @property
        def name(self):
            # returns the name of the predictor
            return self.__name

        @property
        def supportedAlleles(self):
            # returns the supported alleles as strings (without the HLA prefix)
            return self.__alleles

        @property
        def supportedLength(self):
            # returns the supported epitope lengths as iterable
            return self.__supported_length

        @property
        def version(self):
            # returns the version of the predictor
            return self.__version

        # the interface defines a function converting Fred2's HLA allele presentation
        # into an internal presentation used by different methods.
        # for this predictor we won't need it but still have to provide it!
        # the function consumes a list of alleles and converts them into the internally used presentation
        def convert_alleles(self, alleles):
            # we just use the identity function
            return alleles

        # additionally the interface defines a function `predict`
        # that consumes a list of peptides or a single peptide and optionally a list of allele objects
        # this method implements the complete prediction routine
        def predict(self, peptides, alleles=None, binary=False, **kwargs):

            # test whether one peptide or a list
            if isinstance(peptides, basestring):
                peptides = list(peptides)

            # if no alleles are specified do predictions for all supported alleles
            if alleles is None:
                alleles = self.supportedAlleles
            else:
                # filter for supported alleles
                alleles = filter(lambda a: a.name in self.supportedAlleles, alleles)

            result = {}

            # fetch peptides as strings
            peptides = [str(peptide) for peptide in peptides]

            # write peptides temporarily, new line separated
            tmp_input_file = tempfile.NamedTemporaryFile().name
            with open(tmp_input_file, 'wb') as file:
                for peptide in peptides:
                    file.write(peptide + "\n")

            # predict bindings
            for a in alleles:
                result[a] = {}
                tmp_output_file = tempfile.NamedTemporaryFile().name
                mhcnuggets_predict(class_='I',
                                   peptides_path=tmp_input_file,
                                   mhc=a,
                                   output=tmp_output_file + a)

                # read predicted binding affinities back
                with open(tmp_output_file + a, 'rb') as csvfile:
                    reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                    # skip header
                    reader.next()

                    for row in reader:
                        content = row[0].split(',')
                        peptide = content[0]
                        binding_affinity = content[1]
                        if binary:
                            if binding_affinity <= 500:
                                result[a][peptide] = 1.0
                            else:
                                result[a][peptide] = 0.0
                        else:
                            result[a][peptide] = binding_affinity

            # create EpitopePredictionResult object. This is a multi-indexed DataFrame
            # with Peptide and Method as multi-index and alleles as columns
            df_result = EpitopePredictionResult.from_dict(result)
            df_result.index = pandas.MultiIndex.from_tuples([tuple((i, self.name)) for i in df_result.index],
                                                            names=['Seq', 'Method'])
            return df_result
except BadSignatureException:
    print "Class MHCNuggetsPredictor_1 cannot be constructed, because of a bad method signature (predict)"

try:
    class MHCNuggetsPredictor_2(AANNEpitopePrediction):
        """
        Implements MHCNuggets Class II

        .. note::
            Evaluation of machine learning methods to predict peptide binding to MHC Class I proteins
            Rohit Bhattacharya, Ashok Sivakumar, Collin Tokheim, Violeta Beleva Guthrie, Valsamo Anagnostou,
            Victor E. Velculescu, Rachel Karchin (2017) bioRxiv
        """
        __metaclass__ = SignatureCheckerMeta
        __alleles = frozenset(["HLA-DPA1*01:03-DPB1*02:01", "HLA-DPA1*01:03-DPB1*03:01", "HLA-DPA1*01:03-DPB1*04:01",
                               "HLA-DPA1*01:03-DPB1*04:02", "HLA-DPA1*02:01-DPB1*01:01", "HLA-DPA1*02:01-DPB1*05:01",
                               "HLA-DPA1*02:02-DPB1*05:01", "HLA-DPA1*03:01-DPB1*04:02", "HLA-DPB1*01:01",
                               "HLA-DPB1*02:01", "HLA-DPB1*03:01", "HLA-DPB1*04:01", "HLA-DPB1*04:02", "HLA-DPB1*05:01",
                               "HLA-DPB1*09:01", "HLA-DPB1*11:01", "HLA-DPB1*14:01", "HLA-DPB1*20:01", "HLA-DQA1*01:01",
                               "HLA-DQA1*01:01-DQB1*05:01", "HLA-DQA1*01:01-DQB1*05:03", "HLA-DQA1*01:02",
                               "HLA-DQA1*01:02-DQB1*05:01", "HLA-DQA1*01:02-DQB1*05:02", "HLA-DQA1*01:02-DQB1*06:02",
                               "HLA-DQA1*01:02-DQB1*06:04", "HLA-DQA1*01:03-DQB1*03:02", "HLA-DQA1*01:03-DQB1*06:01",
                               "HLA-DQA1*01:03-DQB1*06:03", "HLA-DQA1*01:04-DQB1*05:03", "HLA-DQA1*02:01-DQB1*02:01",
                               "HLA-DQA1*02:01-DQB1*02:02", "HLA-DQA1*02:01-DQB1*03:01", "HLA-DQA1*02:01-DQB1*03:03",
                               "HLA-DQA1*02:01-DQB1*04:02", "HLA-DQA1*03:01", "HLA-DQA1*03:01-DQB1*02:01",
                               "HLA-DQA1*03:01-DQB1*03:01", "HLA-DQA1*03:01-DQB1*03:02", "HLA-DQA1*03:01-DQB1*04:01",
                               "HLA-DQA1*03:02-DQB1*03:01", "HLA-DQA1*03:02-DQB1*03:03", "HLA-DQA1*03:02-DQB1*04:01",
                               "HLA-DQA1*03:03-DQB1*04:02", "HLA-DQA1*04:01-DQB1*04:02", "HLA-DQA1*05:01",
                               "HLA-DQA1*05:01-DQB1*02:01", "HLA-DQA1*05:01-DQB1*03:01", "HLA-DQA1*05:01-DQB1*03:02",
                               "HLA-DQA1*05:01-DQB1*03:03", "HLA-DQA1*05:01-DQB1*04:02", "HLA-DQA1*05:05-DQB1*03:01",
                               "HLA-DQA1*06:01-DQB1*04:02", "HLA-DQB1*02:01", "HLA-DQB1*02:02", "HLA-DQB1*03:01",
                               "HLA-DQB1*03:02", "HLA-DQB1*03:19", "HLA-DQB1*04:02", "HLA-DQB1*05:01", "HLA-DQB1*05:02",
                               "HLA-DQB1*05:03", "HLA-DQB1*06:02", "HLA-DQB1*06:03", "HLA-DQB1*06:04",
                               "HLA-DRA0*10:1-DRB1*01:01", "HLA-DRA0*10:1-DRB1*03:01", "HLA-DRA0*10:1-DRB1*04:01",
                               "HLA-DRA0*10:1-DRB1*04:04", "HLA-DRA0*10:1-DRB1*07:01", "HLA-DRA0*10:1-DRB1*08:01",
                               "HLA-DRA0*10:1-DRB1*09:01", "HLA-DRA0*10:1-DRB1*11:01", "HLA-DRA0*10:1-DRB1*13:01",
                               "HLA-DRA0*10:1-DRB1*14:54", "HLA-DRA0*10:1-DRB1*15:01", "HLA-DRA0*10:1-DRB3*01:01",
                               "HLA-DRA0*10:1-DRB3*02:02", "HLA-DRA0*10:1-DRB3*03:01", "HLA-DRA0*10:1-DRB4*01:03",
                               "HLA-DRA0*10:1-DRB5*01:01", "HLA-DRB1*01:01", "HLA-DRB1*01:02", "HLA-DRB1*01:03",
                               "HLA-DRB1*03:01", "HLA-DRB1*03:02", "HLA-DRB1*03:03", "HLA-DRB1*03:04", "HLA-DRB1*03:05",
                               "HLA-DRB1*04:01", "HLA-DRB1*04:02", "HLA-DRB1*04:03", "HLA-DRB1*04:04", "HLA-DRB1*04:05",
                               "HLA-DRB1*04:06", "HLA-DRB1*04:07", "HLA-DRB1*04:11", "HLA-DRB1*07:01", "HLA-DRB1*08:01",
                               "HLA-DRB1*08:02", "HLA-DRB1*08:03", "HLA-DRB1*08:04", "HLA-DRB1*09:01", "HLA-DRB1*10:01",
                               "HLA-DRB1*11:01", "HLA-DRB1*11:02", "HLA-DRB1*11:03", "HLA-DRB1*11:04", "HLA-DRB1*12:01",
                               "HLA-DRB1*12:02", "HLA-DRB1*13:01", "HLA-DRB1*13:02", "HLA-DRB1*13:03", "HLA-DRB1*13:04",
                               "HLA-DRB1*13:05", "HLA-DRB1*14:01", "HLA-DRB1*14:02", "HLA-DRB1*15:01", "HLA-DRB1*15:02",
                               "HLA-DRB1*15:03", "HLA-DRB1*16:01", "HLA-DRB1*16:02", "HLA-DRB3*01:01", "HLA-DRB3*02:02",
                               "HLA-DRB3*03:01", "HLA-DRB4*01:01", "HLA-DRB4*01:03", "HLA-DRB5*01:01", "HLA-DRB5*01:02"])
        __supported_length = frozenset([5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25])
        __name = "mhcnuggets-class-2"
        __version = "2.0"

        # the interface defines three class properties
        @property
        def name(self):
            # returns the name of the predictor
            return self.__name

        @property
        def supportedAlleles(self):
            # returns the supported alleles as strings (without the HLA prefix)
            return self.__alleles

        @property
        def supportedLength(self):
            # returns the supported epitope lengths as iterable
            return self.__supported_length

        @property
        def version(self):
            # returns the version of the predictor
            return self.__version

        # the interface defines a function converting Fred2's HLA allele presentation
        # into an internal presentation used by different methods.
        # for this predictor we won't need it but still have to provide it!
        # the function consumes a list of alleles and converts them into the internally used presentation
        def convert_alleles(self, alleles):
            return [allele.replace(':', '').replace('*', '') for allele in alleles]

        # additionally the interface defines a function `predict`
        # that consumes a list of peptides or a single peptide and optionally a list
        # of allele objects
        #
        # this method implements the complete prediction routine
        def predict(self, peptides, alleles=None, binary=False, **kwargs):

            # test whether one peptide or a list
            if isinstance(peptides, basestring):
                peptides = list(peptides)

            # if no alleles are specified do predictions for all supported alleles
            if alleles is None:
                alleles = self.supportedAlleles
            else:
                # filter for supported alleles
                alleles = filter(lambda a: a.name in self.supportedAlleles, alleles)

            result = {}

            # fetch peptides as strings
            peptides = [str(peptide) for peptide in peptides]

            # write peptides temporarily, new line separated
            tmp_input_file = tempfile.NamedTemporaryFile().name
            with open(tmp_input_file, 'wb') as file:
                for peptide in peptides:
                    file.write(peptide + "\n")

            # predict bindings
            for a in alleles:
                result[a] = {}
                tmp_output_file = tempfile.NamedTemporaryFile().name

                mhcnuggets_predict(class_='II',
                                   peptides_path=tmp_input_file,
                                   mhc=a,
                                   output=tmp_output_file + a)

                # read predicted binding affinities back
                with open(tmp_output_file + a, 'rb') as csvfile:
                    reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                    # skip header
                    reader.next()

                    for row in reader:
                        content = row[0].split(',')
                        peptide = content[0]
                        binding_affinity = content[1]
                        if binary:
                            if binding_affinity <= 500:
                                result[a][peptide] = 1.0
                            else:
                                result[a][peptide] = 0.0
                        else:
                            result[a][peptide] = binding_affinity

            # create EpitopePredictionResult object. This is a multi-indexed DataFrame
            # with Peptide and Method as multi-index and alleles as columns
            df_result = EpitopePredictionResult.from_dict(result)
            df_result.index = pandas.MultiIndex.from_tuples([tuple((i, self.name)) for i in df_result.index],
                                                            names=['Seq', 'Method'])
            return df_result
except BadSignatureException:
    print "Class MHCNuggetsPredictor_2 cannot be constructed, because of a bad method signature (predict)"

try:
    class MHCFlurryPredictor(AANNEpitopePrediction):
        """
        Implements MHCFlurry

        .. note::
            T. J. O’Donnell, A. Rubinsteyn, M. Bonsack, A. B. Riemer, U. Laserson, and J. Hammerbacher,
             "MHCflurry: Open-Source Class I MHC Binding Affinity Prediction," Cell Systems, 2018.
              Available at: https://www.cell.com/cell-systems/fulltext/S2405-4712(18)30232-1.
        """
        __metaclass__ = SignatureCheckerMeta
        __alleles = frozenset(
            ["HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:02", "HLA-A*02:03", "HLA-A*02:05", "HLA-A*02:06", "HLA-A*02:07",
             "HLA-A*02:11", "HLA-A*02:12", "HLA-A*02:16", "HLA-A*02:17", "HLA-A*02:19", "HLA-A*02:50", "HLA-A*03:01",
             "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02", "HLA-A*24:03", "HLA-A*25:01", "HLA-A*26:01", "HLA-A*26:02",
             "HLA-A*26:03", "HLA-A*29:02", "HLA-A*30:01", 'HLA-A*30:02', "HLA-A*31:01", "HLA-A*32:01", "HLA-A*33:01",
             "HLA-A*66:01", "HLA-A*68:01", "HLA-A*68:02", "HLA-A*68:23", "HLA-A*69:01", "HLA-A*80:01", "HLA-B*07:01",
             "HLA-B*07:02", "HLA-B*08:01", "HLA-B*08:02", "HLA-B*08:03", "HLA-B*14:02", "HLA-B*15:01", "HLA-B*15:02",
             "HLA-B*15:03", "HLA-B*15:09", "HLA-B*15:17", "HLA-B*18:01", "HLA-B*27:02", "HLA-B*27:03", "HLA-B*27:04",
             "HLA-B*27:05", "HLA-B*27:06", "HLA-B*35:01", "HLA-B*35:03", "HLA-B*37:01", "HLA-B*38:01", "HLA-B*39:01",
             "HLA-B*39:06", "HLA-B*40:01", "HLA-B*40:02", "HLA-B*42:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*45:01",
             "HLA-B*46:01", "HLA-B*48:01", "HLA-B*51:01", "HLA-B*53:01", "HLA-B*54:01", "HLA-B*57:01", "HLA-B*58:01",
             "HLA-B*83:01", "HLA-C*03:03", "HLA-C*04:01", "HLA-C*05:01", "HLA-C*06:02", "HLA-C*07:02", "HLA-C*08:02",
             "HLA-C*12:03", "HLA-C*14:02", "HLA-C*15:02"])
        __supported_length = frozenset([8, 9, 10, 11, 12, 13, 14, 15])
        __name = "mhcflurry"
        __version = "1.2.2"

        # the interface defines three class properties
        @property
        def name(self):
            # returns the name of the predictor
            return self.__name

        @property
        def supportedAlleles(self):
            # returns the supported alleles as strings (without the HLA prefix)
            return self.__alleles

        @property
        def supportedLength(self):
            # returns the supported epitope lengths as iterable
            return self.__supported_length

        @property
        def version(self):
            # returns the version of the predictor
            return self.__version

        # the interface defines a function converting Fred2's HLA allele presentation
        # into an internal presentation used by different methods.
        # for this predictor we won't need it but still have to provide it!
        # the function consumes a list of alleles and converts them into the internally used presentation
        def convert_alleles(self, alleles):
            # we just use the identity function
            return alleles

        # additionally the interface defines a function `predict`
        def predict(self, peptides, alleles=None, binary=False, **kwargs):

            # test whether one peptide or a list
            if isinstance(peptides, basestring):
                peptides = list(peptides)

            # if no alleles are specified do predictions for all supported alleles
            if alleles is None:
                alleles = self.supportedAlleles
            else:
                # filter for supported alleles
                alleles = filter(lambda a: a.name in self.supportedAlleles, alleles)

            # test mhcflurry models are available => download if not
            p = subprocess.Popen(['mhcflurry-downloads', 'path', 'models_class1'],
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if p is not 0:
                subprocess.call(['mhcflurry-downloads', 'fetch', 'models_class1'])

            # load model
            predictor = Class1AffinityPredictor.load()

            result = {}

            for a in alleles:
                result[a] = {}
                for p in peptides:
                    seq = p.__str__()
                    binding_affinity = predictor.predict(allele=a, peptides=[seq])
                    if binary:
                        if binding_affinity <= 500:
                            result[a][p] = 1.0
                        else:
                            result[a][p] = 0.0
                    else:
                        result[a][p] = binding_affinity

            # create EpitopePredictionResult object. This is a multi-indexed DataFrame
            # with Peptide and Method as multi-index and alleles as columns
            df_result = EpitopePredictionResult.from_dict(result)
            df_result.index = pandas.MultiIndex.from_tuples([tuple((i, self.name)) for i in df_result.index],
                                                            names=['Seq', 'Method'])
            return df_result
except BadSignatureException:
    print "Class MHCFlurryPredictor cannot be constructed, because of a bad method signature (predict)"
