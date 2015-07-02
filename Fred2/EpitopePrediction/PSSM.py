# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: EpitopePrediction.PSSM
   :synopsis: This module contains all classes for PSSM-based epitope prediction.
.. moduleauthor:: schubert

"""

import itertools
import warnings
import pandas
import math

from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Result import EpitopePredictionResult
from Fred2.Core.Base import AEpitopePrediction


class APSSMEpitopePrediction(AEpitopePrediction):
    """
        Abstract base class for PSSM predictions.

        Implements predict functionality

    """

    def predict(self, peptides, alleles=None, **kwargs):
        """
        Returns predictions for given peptides an alleles. If no alleles are given, predictions for all available models
        are made.

        :param list(Peptide)/Peptide peptides: A single Peptide or a list of Peptides
        :param list(Alleles) alleles: a list of Alleles
        :param kwargs: optional parameter (not used yet)
        :return: Returns a Result object with the prediction results
        """
        def __load_allele_model(allele,length):
            allele_model = "%s_%s_%i"%(self.name, allele, length)
            return getattr( __import__("Fred2.Data.EpitopePSSMMatrices", fromlist=[allele_model]), allele_model)

        if isinstance(peptides, Peptide):
            pep_seqs = {str(peptides):peptides}
        else:
            if any(not isinstance(p, Peptide) for p in peptides):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(p):p for p in peptides}

        if alleles is None:
            al = [Allele("HLA-"+a) for a in self.supportedAlleles]
            allales_string = {conv_a:a for conv_a, a in itertools.izip(self.convert_alleles(al), al)}
        else:
            if isinstance(alleles, Allele):
                alleles = [alleles]
            if any(not isinstance(p, Allele) for p in alleles):
                raise ValueError("Input is not of type Allele")
            allales_string ={conv_a:a for conv_a, a in itertools.izip(self.convert_alleles(alleles),alleles)}

        #group peptides by length and
        result = {}
        for length, peps in itertools.groupby(pep_seqs.iterkeys(), key= lambda x: len(x)):
            peps = list(peps)
            #dynamicaly import prediction PSSMS for alleles and predict
            if length not in self.supportedLength:
                warnings.warn("Peptide length of %i is not supported by %s"%(length,self.name))
                continue

            for a in allales_string.keys():
                try:
                    pssm = __load_allele_model(a, length)
                except AttributeError:
                    warnings.warn("No model found for %s with length %i"%(allales_string[a], length))
                    continue

                result[allales_string[a]] = {}
                ##here is the prediction and result object missing##
                for p in peps:
                    score = sum(pssm[i].get(p[i], 0.0) for i in xrange(length))+pssm.get(-1,{}).get("con", 0)
                    result[allales_string[a]][pep_seqs[p]] = score
                    #print a, score, result

        if not result:
            raise ValueError("No predictions could be made with " +self.name+" for given input. Check your"
                             "epitope length and HLA allele combination.")


        df_result = EpitopePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i,self.name)) for i in df_result.index],
                                                        names=['Seq','Method'])
        return df_result


class Syfpeithi(APSSMEpitopePrediction):
    """
        Represents the Syfpeithi PSSM predictor
    """
    __alleles = frozenset(['B*15:10', 'B*41:01', 'B*37:01', 'B*27:05', 'B*38:01', 'A*02:01', 'B*47:01', 'A*26:01', 'B*37:01',
                 'DRB1*11:01', 'B*50:01', 'B*07:02', 'A*68:01', 'A*24:02', 'DRB1*15:01', 'B*15:01', 'B*45:01',
                 'A*11:01', 'A*03:01', 'B*40:01', 'DRB1*03:01', 'B*39:01', 'DRB1*01:01', 'B*51:01', 'B*39:02',
                 'B*08:01', 'B*18:01', 'B*44:02', 'B*49:01', 'DRB1*07:01', 'B*14:02', 'A*01:01'])
    __supported_length = frozenset([8, 9, 10, 11])
    __name = "syfpeithi"

    @property
    def name(self):
        return self.__name

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__supported_length

    def convert_alleles(self, alleles):
        return ["%s_%s%s"%(a.locus, a.supertype, a.subtype) for a in alleles]

    def predict(self, peptides, alleles=None, **kwargs):
        #with this implementation customizations of prediction algorithm is still possible
        return super(Syfpeithi, self).predict(peptides, alleles=alleles, **kwargs)


class BIMAS(APSSMEpitopePrediction):
    """
        Represents the BIMAS PSSM predictor
    """

    __alleles = frozenset(['B*04:01', 'A*31:01', 'B*58:01', 'C*06:02', 'A*03:01', 'B*35:01',
                 'B*35:01', 'B*15:01', 'A*02:05', 'B*27:05', 'B*27:05', 'A*33:02',
                 'B*39:01', 'B*38:01', 'B*40:', 'A*24:02', 'B*51:01', 'B*07:02', 'B*08:01',
                 'B*51:02', 'B*40:06', 'B*40:06', 'B*51:02', 'B*37:01', 'A*11:01',
                 'B*08:01', 'B*44:03', 'A*68:01', 'B*51:03', 'B*52:01', 'A*02:01',
                 'A*01:01', 'C*07:02', 'C*03:01', 'B*40:01', 'B*51:01', 'B*39:02',
                 'B*52:01', 'C*04:01', 'B*27:02', 'B*39:01'])
    __supported_length = frozenset([8, 9])
    __name = "bimas"

    @property
    def name(self):
        return self.__name

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__supported_length

    def convert_alleles(self, alleles):
        return ["%s_%s%s"%(a.locus, a.supertype, a.subtype) for a in alleles]

    def predict(self, peptides, alleles=None, **kwargs):
        #with this implementation customizations of prediction algorithm is still possible
        return super(BIMAS, self).predict(peptides, alleles=alleles, **kwargs).applymap(lambda x: math.pow(math.e,x))


class Epidemix(APSSMEpitopePrediction):
    __alleles = frozenset(['B*27', 'A*11:01', 'B*27:05', 'B*07', 'B*27', 'A*01', 'B*44', 'A*03',
                 'A*25', 'B*37:01', 'A*02:01', 'A*02:01', 'B*18:01', 'B*18:01', 'A*03',
                 'A*24', 'A*25', 'A*02:01', 'A*11:01', 'A*24:02', 'B*08', 'B*08',
                 'B*51:01', 'B*51:01'])
    __supported_length = frozenset([9, 10, 8, 11])
    __name = "epidemix"


    @property
    def name(self):
        return self.__name

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__supported_length

    def convert_alleles(self, alleles):
        return ["%s_%s%s"%(a.locus, a.supertype, a.subtype) for a in alleles]

    def predict(self, peptides, alleles=None, **kwargs):
        #with this implementation customizations of prediction algorithm is still possible
        return super(Epidemix, self).predict(peptides, alleles=alleles, **kwargs)


class Hammer(APSSMEpitopePrediction):
    __alleles = frozenset(['DRB1*07:03', 'DRB1*07:01', 'DRB1*11:28', 'DRB1*11:21', 'DRB1*11:20', 'DRB1*04:26', 'DRB1*04:23',
                 'DRB1*04:21', 'DRB5*01:05:', 'DRB1*08:17', 'DRB1*13:05', 'DRB1*13:04', 'DRB1*13:07', 'DRB1*13:01',
                 'DRB1*13:02', 'DRB1*08:04', 'DRB1*08:06', 'DRB1*08:01', 'DRB1*01:01', 'DRB1*01:02', 'DRB1*08:02',
                 'DRB1*13:11', 'DRB1*03:11', 'DRB1*11:07', 'DRB1*11:06', 'DRB1*11:04', 'DRB1*11:02', 'DRB1*11:01',
                 'DRB1*04:08', 'DRB1*04:01', 'DRB1*04:02', 'DRB1*04:05', 'DRB1*04:04', 'DRB1*13:23', 'DRB1*13:22',
                 'DRB1*13:21', 'DRB1*13:27', 'DRB1*08:13', 'DRB1*13:28', 'DRB1*03:06', 'DRB1*03:07', 'DRB1*03:05',
                 'DRB1*11:14', 'DRB1*03:01', 'DRB1*15:02', 'DRB1*15:01', 'DRB1*15:06', 'DRB1*03:08', 'DRB1*03:09',
                 'DRB1*04:10', 'DRB5*01:01'])
    __supported_length = frozenset([9])
    __name = "hammer"


    @property
    def name(self):
        return self.__name

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__supported_length

    def convert_alleles(self, alleles):
        return ["%s_%s%s"%(a.locus, a.supertype, a.subtype) for a in alleles]

    def predict(self, peptides, alleles=None, **kwargs):
        #with this implementation customizations of prediction algorithm is still possible
        return super(Hammer, self).predict(peptides, alleles=alleles, **kwargs)


class SMM(APSSMEpitopePrediction):
    """
    Implements IEDBs SMM PSSM method
    """

    __alleles = frozenset(['B*27:20', 'B*83:01', 'A*32:15', 'B*15:17', 'B*40:13', 'A*24:02', 'A*24:03', 'B*53:01', 'B*15:01',
                 'B*27:05', 'B*42:01', 'B*39:01', 'B*38:01', 'A*23:01', 'A*25:01', 'C*04:01', 'A*29:02', 'A*02:06',
                 'A*02:01', 'A*02:02', 'A*02:03', 'A*26:02', 'A*26:03', 'A*26:01', 'C*03:03', 'E*01:01', 'E*01:03',
                 'B*58:01', 'A*31:01', 'C*06:02', 'B*07:02', 'A*66:01', 'B*57:01', 'A*68:01', 'A*68:02', 'C*14:02',
                 'B*35:01', 'B*15:09', 'B*35:03', 'A*80:01', 'B*15:03', 'B*15:02', 'A*32:01', 'A*02:50', 'A*32:07',
                 'B*58:02', 'A*69:01', 'A*68:23', 'A*11:01', 'A*03:01', 'B*73:01', 'B*40:01', 'B*44:03', 'B*46:01',
                 'B*40:02', 'C*12:03', 'B*44:02', 'A*30:01', 'A*02:19', 'A*30:02', 'A*02:17', 'A*02:16', 'B*51:01',
                 'B*45:01', 'A*02:12', 'A*02:11', 'B*54:01', 'B*08:01', 'B*18:01', 'B*08:03', 'B*08:02', 'C*05:01',
                 'C*15:02', 'A*33:01', 'B*14:02', 'C*07:01', 'B*48:01', 'B*15:42', 'C*07:02', 'A*01:01', 'C*08:02'])
    __supported_length = frozenset([8, 9, 10, 11])
    __name = "smm"
    @property
    def name(self):
        return self.__name

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__supported_length

    def convert_alleles(self, alleles):
        return ["%s_%s_%s"%(a.locus, a.supertype, a.subtype) for a in alleles]

    def predict(self, peptides, alleles=None, **kwargs):
        #with this implementation customizations of prediction algorithm is still possible
        #In IEDB scripts score is taken to the base 10**score
        return super(SMM, self).predict(peptides, alleles=alleles, **kwargs).applymap(lambda x: math.pow(10,x))


class SMMPMBEC(APSSMEpitopePrediction):
    """
    Implements IEDBs SMMPMBEC PSSM method
    """

    __alleles = frozenset(['A*01:01', 'A*02:01', 'A*02:02', 'A*02:03', 'A*02:06', 'A*02:11', 'A*02:12', 'A*02:16', 'A*02:17',
                 'A*02:19', 'A*02:50', 'A*03:01', 'A*11:01', 'A*23:01', 'A*24:02', 'A*24:03', 'A*25:01', 'A*26:01',
                 'A*26:02', 'A*26:03', 'A*29:02', 'A*30:01', 'A*30:02', 'A*31:01', 'A*32:01', 'A*32:07', 'A*32:15',
                 'A*33:01', 'A*66:01', 'A*68:01', 'A*68:02', 'A*68:23', 'A*69:01', 'A*80:01', 'B*07:02', 'B*08:01',
                 'B*08:02', 'B*08:03', 'B*14:02', 'B*15:01', 'B*15:02', 'B*15:03', 'B*15:09', 'B*15:17', 'B*15:42',
                 'B*18:01', 'B*27:03', 'B*27:05', 'B*27:20', 'B*35:01', 'B*35:03', 'B*38:01', 'B*39:01', 'B*40:01',
                 'B*40:02', 'B*40:13', 'B*42:01', 'B*44:02', 'B*44:03', 'B*45:01', 'B*45:06', 'B*46:01', 'B*48:01',
                 'B*51:01', 'B*53:01', 'B*54:01', 'B*57:01', 'B*58:01', 'B*58:02', 'B*73:01', 'B*83:01', 'C*03:03',
                 'C*04:01', 'C*05:01', 'C*06:02', 'C*07:01', 'C*07:02', 'C*08:02', 'C*12:03', 'C*14:02', 'C*15:02',
                 'E*01:01', 'E*01:03'])
    __supported_length = frozenset([8, 9, 10, 11])
    __name = "smmpmbec"

    @property
    def name(self):
        return self.__name

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__supported_length

    def convert_alleles(self, alleles):
        return ["%s_%s_%s"%(a.locus, a.supertype, a.subtype) for a in alleles]

    def predict(self, peptides, alleles=None, **kwargs):
        #with this implementation customizations of prediction algorithm is still possible
        #In IEDB scripts score is taken to the base 10**score
        return super(SMMPMBEC, self).predict(peptides, alleles=alleles, **kwargs).applymap(lambda x: math.pow(10,x))


class ARB(APSSMEpitopePrediction):
    """
    Implements IEDBs ARB method
    """

    __alleles = frozenset(['A*01:01', 'A*02:01', 'A*02:02', 'A*02:03', 'A*02:06', 'A*02:11', 'A*02:12', 'A*02:16', 'A*02:19',
                 'A*02:50', 'A*03:01', 'A*11:01', 'A*23:01', 'A*24:02', 'A*24:03', 'A*25:01', 'A*26:01', 'A*26:02',
                 'A*26:03', 'A*29:02', 'A*30:01', 'A*30:02', 'A*31:01', 'A*32:01', 'A*33:01', 'A*68:01', 'A*68:02',
                 'A*69:01', 'A*80:01', 'B*07:02', 'B*08:01', 'B*08:02', 'B*08:03', 'B*15:01', 'B*15:02', 'B*15:03',
                 'B*15:09', 'B*15:17', 'B*18:01', 'B*27:03', 'B*27:05', 'B*35:01', 'B*38:01', 'B*39:01', 'B*40:01',
                 'B*40:02', 'B*44:02', 'B*44:03', 'B*45:01', 'B*46:01', 'B*48:01', 'B*51:01', 'B*53:01', 'B*54:01',
                 'B*57:01', 'B*58:01', 'B*73:01'])
    __supported_length = frozenset([8, 9, 10, 11])
    __name = "arb"

    @property
    def name(self):
        return self.__name

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__supported_length

    def convert_alleles(self, alleles):
        return ["%s_%s%s"%(a.locus, a.supertype, a.subtype) for a in alleles]

    def predict(self, peptides, alleles=None, **kwargs):
        #with this implementation customizations of prediction algorithm is still possible
        #In IEDB scripts score is taken to the base 10**score
        #TODO: IEDB calculates ARB score as follows.... needs some restructuring to be able to do this
        #    score/=-self.length
        #    score-=self.intercept
        #    score/=self.slope
        #    score=math.pow(10,score)
        #    if score < 0.0001:    # Cap predictable values
        #        score = 0.0001
        #    elif score > 1e6:
        #        score = 1e6
        return super(ARB, self).predict(peptides, alleles=alleles, **kwargs).applymap(lambda x: math.pow(10, x))


class ComblibSidney2008(APSSMEpitopePrediction):
    """
    Implements IEDBs Comblib_Sidney2008 PSSM method
    """

    __alleles = frozenset(['B*35:01', 'B*51:01', 'B*54:01', 'B*58:02', 'A*02:01', 'A*68:02', 'B*27:05', 'B*08:01', 'B*07:02',
                 'A*32:01', 'B*53:01', 'A*30:01', 'B*15:03', 'B*15:01', 'B*58:01'])
    __supported_length = frozenset([9])
    __name = "comblibSidney"

    @property
    def name(self):
        return self.__name

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__supported_length

    def convert_alleles(self, alleles):
        return ["%s_%s%s"%(a.locus, a.supertype, a.subtype) for a in alleles]

    def predict(self, peptides, alleles=None, **kwargs):
        #with this implementation customizations of prediction algorithm is still possible
        #In IEDB scripts score is taken to the base 10**score
        return super(ComblibSidney2008, self).predict(peptides,
                                                      alleles=alleles, **kwargs).applymap(lambda x: math.pow(10, x))


class TEPITOPEpan(APSSMEpitopePrediction):
    """
    Implements TEPITOPEpan
    TEPITOPEpan: Extending TEPITOPE for Peptide Binding Prediction Covering over 700 HLA-DR Molecules
                Lianming Zhang , Yiqing Chen , Hau-San Wong, Shuigeng Zhou, Hiroshi Mamitsuka, Shanfeng Zhu
    """

    __alleles = frozenset(['DRB1*01:01', 'DRB1*01:02', 'DRB1*01:03', 'DRB1*01:04', 'DRB1*01:05', 'DRB1*01:06', 'DRB1*01:07',
                 'DRB1*01:08', 'DRB1*01:09', 'DRB1*01:10', 'DRB1*01:11', 'DRB1*01:12', 'DRB1*01:13', 'DRB1*01:14',
                 'DRB1*01:15', 'DRB1*01:16', 'DRB1*01:17', 'DRB1*01:18', 'DRB1*01:19', 'DRB1*01:20', 'DRB1*01:21',
                 'DRB1*01:22', 'DRB1*01:23', 'DRB1*01:24', 'DRB1*01:25', 'DRB1*01:26', 'DRB1*01:27', 'DRB1*01:28',
                 'DRB1*01:29', 'DRB1*01:30', 'DRB1*01:31', 'DRB1*01:32', 'DRB1*01:34', 'DRB1*01:35', 'DRB1*01:36',
                 'DRB1*03:01', 'DRB1*03:02', 'DRB1*03:03', 'DRB1*03:04', 'DRB1*03:05', 'DRB1*03:06', 'DRB1*03:07',
                 'DRB1*03:08', 'DRB1*03:09', 'DRB1*03:10', 'DRB1*03:11', 'DRB1*03:12', 'DRB1*03:13', 'DRB1*03:14',
                 'DRB1*03:15', 'DRB1*03:16', 'DRB1*03:17', 'DRB1*03:18', 'DRB1*03:19', 'DRB1*03:20', 'DRB1*03:21',
                 'DRB1*03:22', 'DRB1*03:23', 'DRB1*03:24', 'DRB1*03:25', 'DRB1*03:26', 'DRB1*03:27', 'DRB1*03:28',
                 'DRB1*03:29', 'DRB1*03:30', 'DRB1*03:31', 'DRB1*03:32', 'DRB1*03:33', 'DRB1*03:34', 'DRB1*03:35',
                 'DRB1*03:36', 'DRB1*03:37', 'DRB1*03:38', 'DRB1*03:39', 'DRB1*03:40', 'DRB1*03:41', 'DRB1*03:42',
                 'DRB1*03:43', 'DRB1*03:44', 'DRB1*03:45', 'DRB1*03:46', 'DRB1*03:47', 'DRB1*03:48', 'DRB1*03:49',
                 'DRB1*03:50', 'DRB1*03:51', 'DRB1*03:52', 'DRB1*03:53', 'DRB1*03:54', 'DRB1*03:55', 'DRB1*03:56',
                 'DRB1*03:57', 'DRB1*03:58', 'DRB1*03:59', 'DRB1*03:60', 'DRB1*03:61', 'DRB1*03:62', 'DRB1*03:63',
                 'DRB1*03:64', 'DRB1*04:01', 'DRB1*04:02', 'DRB1*04:03', 'DRB1*04:04', 'DRB1*04:05', 'DRB1*04:06',
                 'DRB1*04:07', 'DRB1*04:08', 'DRB1*04:09', 'DRB1*04:10', 'DRB1*04:11', 'DRB1*04:12', 'DRB1*04:13',
                 'DRB1*04:14', 'DRB1*04:15', 'DRB1*04:16', 'DRB1*04:17', 'DRB1*04:18', 'DRB1*04:19', 'DRB1*04:20',
                 'DRB1*04:21', 'DRB1*04:22', 'DRB1*04:23', 'DRB1*04:24', 'DRB1*04:25', 'DRB1*04:26', 'DRB1*04:27',
                 'DRB1*04:28', 'DRB1*04:29', 'DRB1*04:30', 'DRB1*04:31', 'DRB1*04:32', 'DRB1*04:33', 'DRB1*04:34',
                 'DRB1*04:35', 'DRB1*04:36', 'DRB1*04:37', 'DRB1*04:38', 'DRB1*04:39', 'DRB1*04:40', 'DRB1*04:41',
                 'DRB1*04:42', 'DRB1*04:43', 'DRB1*04:44', 'DRB1*04:45', 'DRB1*04:46', 'DRB1*04:47', 'DRB1*04:48',
                 'DRB1*04:49', 'DRB1*04:50', 'DRB1*04:51', 'DRB1*04:52', 'DRB1*04:53', 'DRB1*04:54', 'DRB1*04:55',
                 'DRB1*04:56', 'DRB1*04:57', 'DRB1*04:58', 'DRB1*04:59', 'DRB1*04:60', 'DRB1*04:61', 'DRB1*04:62',
                 'DRB1*04:63', 'DRB1*04:64', 'DRB1*04:65', 'DRB1*04:66', 'DRB1*04:67', 'DRB1*04:68', 'DRB1*04:69',
                 'DRB1*04:70', 'DRB1*04:71', 'DRB1*04:72', 'DRB1*04:73', 'DRB1*04:74', 'DRB1*04:75', 'DRB1*04:76',
                 'DRB1*04:77', 'DRB1*04:78', 'DRB1*04:79', 'DRB1*04:80', 'DRB1*04:82', 'DRB1*04:83', 'DRB1*04:84',
                 'DRB1*04:85', 'DRB1*04:86', 'DRB1*04:87', 'DRB1*04:88', 'DRB1*04:89', 'DRB1*04:90', 'DRB1*04:91',
                 'DRB1*04:92', 'DRB1*04:93', 'DRB1*04:95', 'DRB1*04:96', 'DRB1*04:97', 'DRB1*04:98', 'DRB1*07:01',
                 'DRB1*07:03', 'DRB1*07:04', 'DRB1*07:05', 'DRB1*07:06', 'DRB1*07:07', 'DRB1*07:08', 'DRB1*07:09',
                 'DRB1*07:11', 'DRB1*07:12', 'DRB1*07:13', 'DRB1*07:14', 'DRB1*07:15', 'DRB1*07:16', 'DRB1*07:17',
                 'DRB1*07:18', 'DRB1*07:19', 'DRB1*07:20', 'DRB1*07:21', 'DRB1*08:01', 'DRB1*08:02', 'DRB1*08:03',
                 'DRB1*08:04', 'DRB1*08:05', 'DRB1*08:06', 'DRB1*08:07', 'DRB1*08:08', 'DRB1*08:09', 'DRB1*08:10',
                 'DRB1*08:11', 'DRB1*08:12', 'DRB1*08:13', 'DRB1*08:14', 'DRB1*08:15', 'DRB1*08:16', 'DRB1*08:17',
                 'DRB1*08:18', 'DRB1*08:19', 'DRB1*08:20', 'DRB1*08:21', 'DRB1*08:22', 'DRB1*08:23', 'DRB1*08:24',
                 'DRB1*08:25', 'DRB1*08:26', 'DRB1*08:27', 'DRB1*08:28', 'DRB1*08:29', 'DRB1*08:30', 'DRB1*08:31',
                 'DRB1*08:32', 'DRB1*08:33', 'DRB1*08:34', 'DRB1*08:35', 'DRB1*08:36', 'DRB1*08:37', 'DRB1*08:38',
                 'DRB1*08:39', 'DRB1*08:40', 'DRB1*08:41', 'DRB1*08:42', 'DRB1*08:43', 'DRB1*08:44', 'DRB1*08:45',
                 'DRB1*09:01', 'DRB1*09:02', 'DRB1*09:03', 'DRB1*09:04', 'DRB1*09:05', 'DRB1*09:06', 'DRB1*09:07',
                 'DRB1*09:08', 'DRB1*09:09', 'DRB1*09:10', 'DRB1*09:11', 'DRB1*09:12', 'DRB1*10:01', 'DRB1*10:02',
                 'DRB1*10:03', 'DRB1*11:01', 'DRB1*11:02', 'DRB1*11:03', 'DRB1*11:04', 'DRB1*11:05', 'DRB1*11:06',
                 'DRB1*11:07', 'DRB1*11:08', 'DRB1*11:09', 'DRB1*11:10', 'DRB1*11:11', 'DRB1*11:12', 'DRB1*11:13',
                 'DRB1*11:14', 'DRB1*11:15', 'DRB1*11:16', 'DRB1*11:17', 'DRB1*11:18', 'DRB1*11:19', 'DRB1*11:20',
                 'DRB1*11:21', 'DRB1*11:22', 'DRB1*11:23', 'DRB1*11:24', 'DRB1*11:25', 'DRB1*11:26', 'DRB1*11:27',
                 'DRB1*11:28', 'DRB1*11:29', 'DRB1*11:30', 'DRB1*11:31', 'DRB1*11:32', 'DRB1*11:33', 'DRB1*11:34',
                 'DRB1*11:35', 'DRB1*11:36', 'DRB1*11:37', 'DRB1*11:38', 'DRB1*11:39', 'DRB1*11:40', 'DRB1*11:41',
                 'DRB1*11:42', 'DRB1*11:43', 'DRB1*11:44', 'DRB1*11:45', 'DRB1*11:46', 'DRB1*11:47', 'DRB1*11:48',
                 'DRB1*11:49', 'DRB1*11:50', 'DRB1*11:51', 'DRB1*11:52', 'DRB1*11:53', 'DRB1*11:54', 'DRB1*11:55',
                 'DRB1*11:56', 'DRB1*11:57', 'DRB1*11:58', 'DRB1*11:59', 'DRB1*11:60', 'DRB1*11:61', 'DRB1*11:62',
                 'DRB1*11:63', 'DRB1*11:64', 'DRB1*11:65', 'DRB1*11:66', 'DRB1*11:67', 'DRB1*11:68', 'DRB1*11:69',
                 'DRB1*11:70', 'DRB1*11:72', 'DRB1*11:73', 'DRB1*11:74', 'DRB1*11:75', 'DRB1*11:76', 'DRB1*11:77',
                 'DRB1*11:78', 'DRB1*11:79', 'DRB1*11:80', 'DRB1*11:81', 'DRB1*11:82', 'DRB1*11:83', 'DRB1*11:84',
                 'DRB1*11:85', 'DRB1*11:86', 'DRB1*11:87', 'DRB1*11:88', 'DRB1*11:89', 'DRB1*11:90', 'DRB1*11:91',
                 'DRB1*11:92', 'DRB1*11:93', 'DRB1*11:94', 'DRB1*11:95', 'DRB1*11:96', 'DRB1*11:97', 'DRB1*11:98',
                 'DRB1*11:99', 'DRB1*12:01', 'DRB1*12:02', 'DRB1*12:03', 'DRB1*12:04', 'DRB1*12:05', 'DRB1*12:06',
                 'DRB1*12:07', 'DRB1*12:08', 'DRB1*12:09', 'DRB1*12:10', 'DRB1*12:11', 'DRB1*12:12', 'DRB1*12:13',
                 'DRB1*12:14', 'DRB1*12:15', 'DRB1*12:16', 'DRB1*12:17', 'DRB1*12:18', 'DRB1*12:19', 'DRB1*12:20',
                 'DRB1*12:21', 'DRB1*12:22', 'DRB1*12:23', 'DRB1*12:25', 'DRB1*12:26', 'DRB1*12:27', 'DRB1*13:01',
                 'DRB1*13:02', 'DRB1*13:03', 'DRB1*13:04', 'DRB1*13:05', 'DRB1*13:06', 'DRB1*13:07', 'DRB1*13:08',
                 'DRB1*13:09', 'DRB1*13:10', 'DRB1*13:11', 'DRB1*13:12', 'DRB1*13:13', 'DRB1*13:14', 'DRB1*13:15',
                 'DRB1*13:16', 'DRB1*13:17', 'DRB1*13:18', 'DRB1*13:19', 'DRB1*13:20', 'DRB1*13:21', 'DRB1*13:22',
                 'DRB1*13:23', 'DRB1*13:24', 'DRB1*13:25', 'DRB1*13:26', 'DRB1*13:27', 'DRB1*13:28', 'DRB1*13:29',
                 'DRB1*13:30', 'DRB1*13:31', 'DRB1*13:32', 'DRB1*13:33', 'DRB1*13:34', 'DRB1*13:35', 'DRB1*13:36',
                 'DRB1*13:37', 'DRB1*13:38', 'DRB1*13:39', 'DRB1*13:40', 'DRB1*13:41', 'DRB1*13:42', 'DRB1*13:43',
                 'DRB1*13:44', 'DRB1*13:45', 'DRB1*13:46', 'DRB1*13:47', 'DRB1*13:48', 'DRB1*13:49', 'DRB1*13:50',
                 'DRB1*13:51', 'DRB1*13:52', 'DRB1*13:53', 'DRB1*13:54', 'DRB1*13:55', 'DRB1*13:56', 'DRB1*13:57',
                 'DRB1*13:58', 'DRB1*13:59', 'DRB1*13:60', 'DRB1*13:61', 'DRB1*13:62', 'DRB1*13:63', 'DRB1*13:64',
                 'DRB1*13:65', 'DRB1*13:66', 'DRB1*13:67', 'DRB1*13:68', 'DRB1*13:69', 'DRB1*13:70', 'DRB1*13:71',
                 'DRB1*13:72', 'DRB1*13:73', 'DRB1*13:74', 'DRB1*13:75', 'DRB1*13:76', 'DRB1*13:77', 'DRB1*13:78',
                 'DRB1*13:79', 'DRB1*13:80', 'DRB1*13:81', 'DRB1*13:82', 'DRB1*13:83', 'DRB1*13:84', 'DRB1*13:85',
                 'DRB1*13:86', 'DRB1*13:87', 'DRB1*13:88', 'DRB1*13:89', 'DRB1*13:90', 'DRB1*13:91', 'DRB1*13:92',
                 'DRB1*13:93', 'DRB1*13:94', 'DRB1*13:95', 'DRB1*13:96', 'DRB1*13:97', 'DRB1*13:98', 'DRB1*13:99',
                 'DRB1*14:01', 'DRB1*14:02', 'DRB1*14:03', 'DRB1*14:04', 'DRB1*14:05', 'DRB1*14:06', 'DRB1*14:07',
                 'DRB1*14:08', 'DRB1*14:09', 'DRB1*14:10', 'DRB1*14:11', 'DRB1*14:12', 'DRB1*14:13', 'DRB1*14:14',
                 'DRB1*14:15', 'DRB1*14:16', 'DRB1*14:17', 'DRB1*14:18', 'DRB1*14:19', 'DRB1*14:20', 'DRB1*14:21',
                 'DRB1*14:22', 'DRB1*14:23', 'DRB1*14:24', 'DRB1*14:25', 'DRB1*14:26', 'DRB1*14:27', 'DRB1*14:28',
                 'DRB1*14:29', 'DRB1*14:30', 'DRB1*14:31', 'DRB1*14:32', 'DRB1*14:33', 'DRB1*14:34', 'DRB1*14:35',
                 'DRB1*14:36', 'DRB1*14:37', 'DRB1*14:38', 'DRB1*14:39', 'DRB1*14:40', 'DRB1*14:41', 'DRB1*14:42',
                 'DRB1*14:43', 'DRB1*14:44', 'DRB1*14:45', 'DRB1*14:46', 'DRB1*14:47', 'DRB1*14:48', 'DRB1*14:49',
                 'DRB1*14:50', 'DRB1*14:51', 'DRB1*14:52', 'DRB1*14:53', 'DRB1*14:54', 'DRB1*14:55', 'DRB1*14:56',
                 'DRB1*14:57', 'DRB1*14:58', 'DRB1*14:59', 'DRB1*14:60', 'DRB1*14:61', 'DRB1*14:62', 'DRB1*14:63',
                 'DRB1*14:64', 'DRB1*14:65', 'DRB1*14:67', 'DRB1*14:68', 'DRB1*14:69', 'DRB1*14:70', 'DRB1*14:71',
                 'DRB1*14:72', 'DRB1*14:73', 'DRB1*14:74', 'DRB1*14:75', 'DRB1*14:76', 'DRB1*14:77', 'DRB1*14:78',
                 'DRB1*14:79', 'DRB1*14:80', 'DRB1*14:81', 'DRB1*14:82', 'DRB1*14:83', 'DRB1*14:84', 'DRB1*14:85',
                 'DRB1*14:86', 'DRB1*14:87', 'DRB1*14:88', 'DRB1*14:89', 'DRB1*14:90', 'DRB1*14:91', 'DRB1*14:93',
                 'DRB1*14:94', 'DRB1*14:95', 'DRB1*14:96', 'DRB1*14:97', 'DRB1*14:98', 'DRB1*14:99', 'DRB1*15:01',
                 'DRB1*15:02', 'DRB1*15:03', 'DRB1*15:04', 'DRB1*15:05', 'DRB1*15:06', 'DRB1*15:07', 'DRB1*15:08',
                 'DRB1*15:09', 'DRB1*15:10', 'DRB1*15:11', 'DRB1*15:12', 'DRB1*15:13', 'DRB1*15:14', 'DRB1*15:15',
                 'DRB1*15:16', 'DRB1*15:18', 'DRB1*15:19', 'DRB1*15:20', 'DRB1*15:21', 'DRB1*15:22', 'DRB1*15:23',
                 'DRB1*15:24', 'DRB1*15:25', 'DRB1*15:26', 'DRB1*15:27', 'DRB1*15:28', 'DRB1*15:29', 'DRB1*15:30',
                 'DRB1*15:31', 'DRB1*15:32', 'DRB1*15:33', 'DRB1*15:34', 'DRB1*15:35', 'DRB1*15:36', 'DRB1*15:37',
                 'DRB1*15:38', 'DRB1*15:39', 'DRB1*15:40', 'DRB1*15:41', 'DRB1*15:42', 'DRB1*15:43', 'DRB1*15:44',
                 'DRB1*15:45', 'DRB1*15:46', 'DRB1*15:47', 'DRB1*15:48', 'DRB1*15:49', 'DRB1*15:51', 'DRB1*15:52',
                 'DRB1*15:53', 'DRB1*15:54', 'DRB1*15:55', 'DRB1*15:56', 'DRB1*15:57', 'DRB1*16:01', 'DRB1*16:02',
                 'DRB1*16:03', 'DRB1*16:04', 'DRB1*16:05', 'DRB1*16:07', 'DRB1*16:08', 'DRB1*16:09', 'DRB1*16:10',
                 'DRB1*16:11', 'DRB1*16:12', 'DRB1*16:14', 'DRB1*16:15', 'DRB1*16:16', 'DRB1*16:17', 'DRB1*16:18',
                 'DRB3*01:01', 'DRB3*01:02', 'DRB3*01:03', 'DRB3*01:04', 'DRB3*01:05', 'DRB3*01:06', 'DRB3*01:07',
                 'DRB3*01:08', 'DRB3*01:09', 'DRB3*01:10', 'DRB3*01:11', 'DRB3*01:12', 'DRB3*01:13', 'DRB3*01:14',
                 'DRB3*01:15', 'DRB3*02:01', 'DRB3*02:02', 'DRB3*02:03', 'DRB3*02:04', 'DRB3*02:05', 'DRB3*02:06',
                 'DRB3*02:07', 'DRB3*02:08', 'DRB3*02:09', 'DRB3*02:10', 'DRB3*02:11', 'DRB3*02:12', 'DRB3*02:13',
                 'DRB3*02:14', 'DRB3*02:15', 'DRB3*02:16', 'DRB3*02:17', 'DRB3*02:18', 'DRB3*02:19', 'DRB3*02:20',
                 'DRB3*02:21', 'DRB3*02:22', 'DRB3*02:23', 'DRB3*02:24', 'DRB3*02:25', 'DRB3*02:26', 'DRB3*02:27',
                 'DRB3*02:28', 'DRB3*03:01', 'DRB3*03:02', 'DRB3*03:03', 'DRB4*01:01', 'DRB4*01:03', 'DRB4*01:04',
                 'DRB4*01:05', 'DRB4*01:06', 'DRB4*01:07', 'DRB4*01:08', 'DRB5*01:01', 'DRB5*01:02', 'DRB5*01:04',
                 'DRB5*01:05', 'DRB5*01:06', 'DRB5*01:07', 'DRB5*01:08', 'DRB5*01:09', 'DRB5*01:11', 'DRB5*01:12',
                 'DRB5*01:13', 'DRB5*01:14', 'DRB5*02:02', 'DRB5*02:03', 'DRB5*02:04', 'DRB5*02:05'])
    __supported_length = frozenset([9])
    __name = "tepitopepan"

    @property
    def name(self):
        return self.__name

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__supported_length

    def convert_alleles(self, alleles):
        return ["%s_%s%s"%(a.locus, a.supertype, a.subtype) for a in alleles]

    def predict(self, peptides, alleles=None, **kwargs):
        return super(TEPITOPEpan, self).predict(peptides, alleles=alleles, **kwargs)
