# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'walzer', 'schubert'

import collections, itertools, warnings,pandas
import math

from Fred2.Core.Allele import Allele
from Fred2.Core.Result import EpitopePredictionResult
from Fred2.Core.Base import AEpitopePrediction


class APSSMEpitopePredictor(AEpitopePrediction):
    """
        Abstract base class for PSSM predictions.

        Implements predict functionality

    """
    def threshold(self, allele):
        return 0.0

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

            #TODO: what if there exists no allele model for this length?
            return getattr( __import__("Fred2.Data.EpitopePSSMMatrices", fromlist=[allele_model]), allele_model)


        if isinstance(peptides, collections.Iterable):
            pep_seqs = {str(p):p for p in peptides}
        else:
            pep_seqs = {str(peptides):peptides}

        if alleles is None:
            al = [Allele("HLA-"+a) for a in self.supportedAlleles]
            allales_string = {conv_a:a for conv_a, a in itertools.izip(self.convert_alleles(al), al)}
        else:
            allales_string ={conv_a:a for conv_a, a in itertools.izip(self.convert_alleles(alleles),alleles)}

        #group peptides by length and
        result = {}
        for length, peps in itertools.groupby(pep_seqs.iterkeys(), key= lambda x: len(x)):
            peps = list(peps)
            #dynamicaly import prediction PSSMS for alleles and predict
            if length not in self.supportedLength:
                warnings.warn("Peptide length of %i not supported"%length, RuntimeWarning)
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
                    score = sum(pssm[i].get(p[i], 0.0) for i in xrange(length))
                    result[allales_string[a]][pep_seqs[p]] = score
                    #print a, score, result

        df_result = EpitopePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i,self.name)) for i in df_result.index],
                                                        names=['Seq','Method'])
        return df_result


class Syfpeithi(APSSMEpitopePredictor):
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


class BIMAS(APSSMEpitopePredictor):
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
        return super(BIMAS, self).predict(peptides, alleles=alleles, **kwargs)


class Epidemix(APSSMEpitopePredictor):
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


class Hammer(APSSMEpitopePredictor):
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


class SMM(APSSMEpitopePredictor):
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
    __offset = {('A*30:02', 10): 4.04038, ('A*02:02', 8): 4.58364, ('A*32:01', 10): 4.33802, ('A*68:02', 8): 4.63836,
                ('A*02:50', 9): 1.93222, ('B*27:05', 8): 4.5795, ('A*23:01', 8): 4.4944, ('A*26:02', 9): 4.20732,
                ('A*29:02', 8): 4.56465, ('B*53:01', 9): 4.62588, ('A*68:01', 9): 4.47649, ('A*24:02', 8): 4.39614,
                ('C*03:03', 9): 2.13061, ('B*44:03', 9): 4.85249, ('C*06:02', 9): 4.48581, ('A*11:01', 9): 5.0053,
                ('A*30:02', 9): 4.09046, ('B*48:01', 9): 5.4796, ('B*08:01', 8): 4.34903, ('A*29:02', 11): 4.15163,
                ('A*03:01', 11): 3.93979, ('B*15:42', 9): 3.70134, ('B*35:03', 9): 5.16319, ('B*35:01', 9): 4.48558,
                ('A*02:06', 9): 4.23955, ('B*44:02', 8): 4.50882, ('A*02:19', 9): 6.31174, ('B*51:01', 11): 4.60685,
                ('B*83:01', 9): 3.87294, ('A*02:03', 10): 4.22725, ('B*39:01', 9): 5.2854, ('B*08:01', 11): 4.36091,
                ('A*31:01', 10): 4.1264, ('B*40:02', 11): 4.46569, ('B*58:02', 9): 4.59052, ('B*35:03', 10): 4.72951,
                ('B*18:01', 10): 4.64684, ('A*02:06', 10): 4.33791, ('B*07:02', 11): 4.59281, ('C*04:01', 9): 4.53642,
                ('B*15:09', 9): 4.84432, ('B*15:02', 9): 3.04223, ('C*07:01', 9): 3.34217, ('A*02:03', 9): 4.65395,
                ('A*02:01', 9): 5.08211, ('A*01:01', 9): 5.07755, ('A*66:01', 9): 4.26663, ('A*26:01', 8): 4.8375,
                ('A*30:01', 10): 4.03394, ('A*68:02', 9): 4.58359, ('A*31:01', 9): 4.59647, ('B*57:01', 9): 5.03382,
                ('A*02:02', 9): 4.16801, ('E*01:01', 9): 4.37622, ('B*27:05', 11): 4.37518, ('B*08:02', 9): 5.2022,
                ('A*23:01', 9): 4.66943, ('B*18:01', 9): 4.75104, ('B*58:01', 11): 4.04406, ('B*53:01', 10): 4.59784,
                ('C*07:02', 9): 2.80617, ('A*02:16', 9): 5.32749, ('A*24:02', 9): 4.56063, ('B*40:01', 8): 4.5203,
                ('C*12:03', 9): 1.48921, ('B*44:03', 11): 4.30511, ('A*11:01', 8): 4.52176, ('A*24:03', 9): 4.43954,
                ('A*30:02', 8): 4.4224, ('B*27:20', 9): 1.18649, ('A*02:02', 10): 4.17451, ('A*33:01', 10): 4.2855,
                ('B*08:01', 9): 5.03255, ('E*01:03', 9): 4.72824, ('B*15:17', 9): 5.37971, ('A*29:02', 10): 4.38028,
                ('A*03:01', 8): 4.58488, ('B*73:01', 9): 4.21773, ('A*24:02', 10): 4.63972, ('B*35:01', 10): 4.6849,
                ('A*02:17', 10): 2.83174, ('B*51:01', 8): 4.60744, ('B*45:01', 10): 4.41769, ('A*11:01', 11): 4.0931,
                ('A*01:01', 10): 4.7508, ('C*14:02', 9): 2.29562, ('A*32:15', 9): 2.16325, ('B*08:03', 9): 4.6631,
                ('C*08:02', 9): 3.50216, ('B*40:02', 10): 4.66289, ('A*25:01', 9): 6.38608, ('A*02:06', 11): 4.05032,
                ('B*07:02', 8): 4.66201, ('B*54:01', 10): 4.36032, ('B*45:01', 9): 4.55821, ('C*05:01', 9): 3.75267,
                ('A*02:03', 8): 4.37112, ('A*02:01', 8): 4.31021, ('A*80:01', 9): 5.91501, ('B*44:03', 8): 4.45534,
                ('A*26:01', 11): 4.32539, ('A*68:02', 11): 4.35192, ('B*57:01', 10): 4.36541, ('B*40:02', 9): 4.69586,
                ('B*15:03', 9): 3.79609, ('A*02:17', 9): 2.99966, ('A*69:01', 9): 4.95386, ('A*02:12', 9): 5.39669,
                ('B*27:05', 10): 4.50327, ('A*02:01', 11): 4.10891, ('A*23:01', 10): 4.48671, ('B*18:01', 8): 4.36736,
                ('B*58:01', 10): 4.20625, ('B*53:01', 11): 4.37161, ('B*40:01', 9): 5.1879, ('B*15:01', 9): 4.76,
                ('C*15:02', 9): 3.65542, ('A*30:02', 11): 3.8759, ('A*02:11', 9): 4.66433, ('A*02:02', 11): 3.8944,
                ('A*33:01', 9): 4.42505, ('B*14:02', 9): 4.30605, ('B*08:01', 10): 4.32962, ('B*27:05', 9): 4.87457,
                ('B*44:03', 10): 4.83507, ('B*18:01', 11): 4.33825, ('A*29:02', 9): 4.45576, ('A*03:01', 9): 5.2158,
                ('B*53:01', 8): 4.61839, ('A*68:01', 10): 4.79599, ('B*54:01', 8): 4.64673, ('A*24:02', 11): 4.13708,
                ('B*35:01', 11): 4.19065, ('B*44:02', 10): 4.55889, ('A*68:02', 10): 4.62831, ('B*51:01', 9): 4.98563,
                ('B*15:01', 10): 3.66225, ('A*11:01', 10): 5.16585, ('A*03:01', 10): 4.7022, ('B*35:01', 8): 4.72662,
                ('A*02:06', 8): 4.36214, ('B*07:02', 9): 5.46489, ('B*44:02', 9): 5.1306, ('B*54:01', 9): 4.49688,
                ('B*51:01', 10): 4.89605, ('B*45:01', 8): 4.60845, ('A*68:23', 9): 1.43775, ('A*02:03', 11): 3.84418,
                ('A*01:01', 11): 4.63802, ('A*26:01', 10): 4.84669, ('B*42:01', 9): 3.85875, ('A*31:01', 11): 4.10452,
                ('B*57:01', 11): 3.8374, ('B*40:02', 8): 4.46963, ('B*15:03', 10): 3.50047, ('B*46:01', 9): 6.47228,
                ('B*07:02', 10): 4.55182, ('B*38:01', 9): 5.63807, ('A*32:07', 9): 1.37905, ('A*32:01', 9): 5.04594,
                ('A*02:01', 10): 4.63048, ('B*40:13', 9): 1.82293, ('A*23:01', 11): 4.16437, ('A*01:01', 8): 4.51394,
                ('B*58:01', 9): 5.0683, ('A*26:01', 9): 5.28461, ('A*30:01', 9): 4.25797, ('A*26:03', 9): 4.73079,
                ('B*42:01', 10): 4.05963, ('B*40:01', 10): 4.75803, ('B*45:01', 11): 4.33229}
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
        r_df = super(SMM, self).predict(peptides, alleles=alleles, **kwargs)
        for i in r_df.index:
            for c in r_df.columns:
                r_df.loc[i,c] = math.pow(10,  r_df.loc[i,c]+self.__offset[c.name,len(i[0])])
        return r_df


class SMMPMBEC(APSSMEpitopePredictor):
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
    __offset = {('A*30:02', 10): 4.11588, ('A*02:02', 8): 4.58388, ('A*32:01', 10): 4.53995, ('A*68:02', 8): 4.73857,
                ('A*02:50', 9): 2.12727, ('B*27:05', 8): 4.73607, ('A*23:01', 8): 4.37966, ('A*26:02', 9): 4.15094,
                ('A*29:02', 8): 4.54295, ('B*53:01', 9): 4.52396, ('A*68:01', 9): 4.50885, ('A*24:02', 8): 4.43431,
                ('C*03:03', 9): 2.40677, ('B*44:03', 9): 4.74672, ('C*06:02', 9): 4.45257, ('A*11:01', 9): 4.98403,
                ('A*30:02', 9): 4.08695, ('B*48:01', 9): 5.11227, ('B*08:01', 8): 4.08066, ('A*29:02', 11): 4.36087,
                ('A*03:01', 11): 4.01267, ('B*15:42', 9): 3.73087, ('B*35:03', 9): 5.21287, ('B*35:01', 9): 4.46264,
                ('A*02:06', 9): 4.27797, ('B*44:02', 8): 4.60445, ('A*02:19', 9): 6.00028, ('B*51:01', 11): 4.43461,
                ('B*83:01', 9): 3.92461, ('A*02:03', 10): 4.23339, ('B*39:01', 9): 5.21899, ('B*08:01', 11): 4.2866,
                ('A*31:01', 10): 4.18421, ('B*40:02', 11): 4.17244, ('B*58:02', 9): 4.3782, ('B*35:03', 10): 4.67913,
                ('B*18:01', 10): 4.65015, ('A*02:06', 10): 4.43229, ('B*07:02', 11): 4.33356, ('C*04:01', 9): 4.55129,
                ('B*15:09', 9): 4.56203, ('B*15:02', 9): 3.4263, ('C*07:01', 9): 3.29367, ('A*02:03', 9): 4.68067,
                ('A*02:01', 9): 5.08043, ('A*01:01', 9): 5.08242, ('A*66:01', 9): 3.97586, ('A*26:01', 8): 4.77635,
                ('A*30:01', 10): 3.96982, ('A*68:02', 9): 4.59764, ('A*31:01', 9): 4.50992, ('B*57:01', 9): 4.87661,
                ('A*02:02', 9): 4.20228, ('E*01:01', 9): 4.25938, ('B*27:05', 11): 4.37011, ('B*08:02', 9): 5.27892,
                ('A*23:01', 9): 4.65717, ('B*18:01', 9): 4.76856, ('B*58:01', 11): 4.09535, ('B*53:01', 10): 4.5361,
                ('C*07:02', 9): 2.88663, ('A*02:16', 9): 5.30329, ('A*24:02', 9): 4.56396, ('B*40:01', 8): 4.58328,
                ('C*12:03', 9): 1.55043, ('B*44:03', 11): 4.27628, ('A*11:01', 8): 4.47965, ('A*24:03', 9): 4.37328,
                ('A*30:02', 8): 4.4448, ('B*27:20', 9): 1.48487, ('A*02:02', 10): 4.16322, ('A*33:01', 10): 4.36026,
                ('B*08:01', 9): 4.94213, ('E*01:03', 9): 4.61031, ('B*15:17', 9): 5.0693, ('A*29:02', 10): 4.50653,
                ('A*03:01', 8): 4.58093, ('B*73:01', 9): 4.36786, ('A*24:02', 10): 4.72765, ('B*35:01', 10): 4.75615,
                ('A*02:17', 10): 2.90381, ('B*51:01', 8): 4.58727, ('B*45:01', 10): 4.39293, ('A*11:01', 11): 4.0985,
                ('A*01:01', 10): 4.75854, ('C*14:02', 9): 2.28569, ('A*32:15', 9): 2.14452, ('B*08:03', 9): 4.84319,
                ('C*08:02', 9): 3.58102, ('B*40:02', 10): 4.52283, ('A*25:01', 9): 5.45384, ('A*02:06', 11): 4.12927,
                ('B*07:02', 8): 4.64764, ('B*54:01', 10): 4.3502, ('B*45:01', 9): 4.53962, ('C*05:01', 9): 3.57027,
                ('A*02:03', 8): 4.2909, ('A*02:01', 8): 4.26183, ('A*80:01', 9): 5.65756, ('B*44:03', 8): 4.52063,
                ('A*26:01', 11): 4.34928, ('A*68:02', 11): 4.32633, ('B*57:01', 10): 4.36574, ('B*40:02', 9): 4.73887,
                ('B*15:03', 9): 3.886, ('A*02:17', 9): 3.43858, ('A*69:01', 9): 4.90648, ('A*02:12', 9): 5.28683,
                ('B*27:05', 10): 4.5218, ('A*02:01', 11): 4.25318, ('A*23:01', 10): 4.47986, ('B*18:01', 8): 4.38825,
                ('B*58:01', 10): 4.22328, ('B*53:01', 11): 4.30256, ('B*40:01', 9): 5.20413, ('B*15:01', 9): 4.72595,
                ('C*15:02', 9): 3.65181, ('A*30:02', 11): 3.86981, ('A*02:11', 9): 4.72423, ('A*02:02', 11): 3.96545,
                ('A*33:01', 9): 4.4772, ('B*14:02', 9): 4.24509, ('B*08:01', 10): 4.33252, ('B*27:05', 9): 4.85457,
                ('B*44:03', 10): 4.87312, ('B*18:01', 11): 4.30021, ('A*29:02', 9): 4.46097, ('A*03:01', 9): 5.19736,
                ('B*53:01', 8): 4.585, ('A*68:01', 10): 4.73869, ('B*54:01', 8): 4.56204, ('A*24:02', 11): 4.15411,
                ('B*35:01', 11): 4.27325, ('B*44:02', 10): 4.54079, ('A*68:02', 10): 4.53515, ('B*51:01', 9): 4.95522,
                ('B*15:01', 10): 3.80964, ('A*11:01', 10): 5.01083, ('A*03:01', 10): 4.73507, ('B*35:01', 8): 4.80159,
                ('A*02:06', 8): 4.38399, ('B*07:02', 9): 5.45316, ('B*44:02', 9): 4.85468, ('B*54:01', 9): 4.48133,
                ('B*51:01', 10): 4.89945, ('B*45:01', 8): 4.58524, ('A*68:23', 9): 1.64682, ('A*02:03', 11): 4.02827,
                ('A*01:01', 11): 4.3952, ('A*26:01', 10): 4.94037, ('B*42:01', 9): 3.88483, ('A*31:01', 11): 4.07222,
                ('B*57:01', 11): 3.71407, ('B*40:02', 8): 4.46747, ('B*15:03', 10): 3.62555, ('B*46:01', 9): 5.83693,
                ('B*07:02', 10): 4.60063, ('B*38:01', 9): 5.51229, ('A*32:07', 9): 1.46981, ('A*32:01', 9): 4.95196,
                ('A*02:01', 10): 4.58782, ('B*40:13', 9): 1.95721, ('A*23:01', 11): 4.10391, ('A*01:01', 8): 4.52761,
                ('B*58:01', 9): 4.98059, ('A*26:01', 9): 5.30881, ('A*30:01', 9): 4.23841, ('A*26:03', 9): 4.43222,
                ('B*42:01', 10): 4.04951, ('B*40:01', 10): 4.74752, ('B*45:01', 11): 4.38364}
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
        r_df = super(SMMPMBEC, self).predict(peptides, alleles=alleles, **kwargs)
        for i in r_df.index:
            for c in r_df.columns:
                r_df.loc[i,c] = math.pow(10,  r_df.loc[i,c]+self.__offset[c.name,len(i[0])])
        return r_df


class ARB(APSSMEpitopePredictor):
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


class ComblibSidney2008(APSSMEpitopePredictor):
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


class TEPITOPEpan(APSSMEpitopePredictor):
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