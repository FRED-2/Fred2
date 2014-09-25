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
            allales_string ={conv_a:a.name for conv_a, a in itertools.izip(self.convert_alleles(alleles),alleles)}

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
                except ImportError:
                    warnings.warn("No model found for %s with length %i"%(allales_string[a], length))
                    continue

                result[allales_string[a]] = {}
                ##here is the prediction and result object missing##
                for p in peps:
                    score = sum(pssm[i][p[i]] for i in xrange(length))
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
    __alleles = ['B*15:10', 'B*41:01', 'B*37:01', 'B*27:05', 'B*38:01', 'A*02:01', 'B*47:01', 'A*26:01', 'B*3701',
                 'DRB1*11:01', 'B*50:01', 'B*07:02', 'A*68:01', 'A*24:02', 'DRB1*15:01', 'B*15:01', 'B*45:01',
                 'A*11:01', 'A*03:01', 'B*40:01', 'DRB1*03:01', 'B*39:01', 'DRB1*01:01', 'B*51:01', 'B*39:02',
                 'B*08:01', 'B*18:01', 'B*44:02', 'B*49:01', 'DRB1*07:01', 'B*14:02', 'A*01:01']
    __supported_length = [8, 9, 10, 11]
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

    __alleles = ['B*04:01', 'A*31:01', 'B*58:01', 'Cw*06:02', 'A*03:01', 'B*35:01',
                 'B*35:01', 'B*15:01', 'A*02:05', 'B*27:05', 'B*27:05', 'A*33:02',
                 'B*39:01', 'B*38:01', 'B*40:', 'A*24:02', 'B*51:01', 'B*07:02', 'B*08:01',
                 'B*51:02', 'B*40:06', 'B*40:06', 'B*51:02', 'B*37:01', 'A*11:01',
                 'B*08:01', 'B*44:03', 'A*68:01', 'B*51:03', 'B*52:01', 'A*02:01',
                 'A*01:01', 'C*07:02', 'C*03:01', 'B*40:01', 'B*51:01', 'B*39:02',
                 'B*52:01', 'C*04:01', 'B*27:02', 'B*39:01']
    __supported_length = [8, 9]
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
    __alleles = ['B*27', 'A*11:01', 'B*27:05', 'B*07', 'B*27', 'A*01', 'B*44', 'A*03',
                 'A*25', 'B*37:01', 'A*02:01', 'A*02:01', 'B*18:01', 'B*18:01', 'A*03',
                 'A*24', 'A*25', 'A*02:01', 'A*11:01', 'A*24:02', 'B*08', 'B*08',
                 'B*51:01', 'B*51:01']
    __supported_length = [9, 10, 8, 11]
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
    __alleles = ['DRB1*07:03', 'DRB1*07:01', 'DRB1*11:28', 'DRB1*11:21', 'DRB1*11:20', 'DRB1*04:26', 'DRB1*04:23',
                 'DRB1*04:21', 'DRB5*01:05:', 'DRB1*08:17', 'DRB1*13:05', 'DRB1*13:04', 'DRB1*13:07', 'DRB1*13:01',
                 'DRB1*13:02', 'DRB1*08:04', 'DRB1*08:06', 'DRB1*08:01', 'DRB1*01:01', 'DRB1*01:02', 'DRB1*08:02',
                 'DRB1*13:11', 'DRB1*03:11', 'DRB1*11:07', 'DRB1*11:06', 'DRB1*11:04', 'DRB1*11:02', 'DRB1*11:01',
                 'DRB1*04:08', 'DRB1*04:01', 'DRB1*04:02', 'DRB1*04:05', 'DRB1*04:04', 'DRB1*13:23', 'DRB1*13:22',
                 'DRB1*13:21', 'DRB1*13:27', 'DRB1*08:13', 'DRB1*13:28', 'DRB1*03:06', 'DRB1*03:07', 'DRB1*03:05',
                 'DRB1*11:14', 'DRB1*03:01', 'DRB1*15:02', 'DRB1*15:01', 'DRB1*15:06', 'DRB1*03:08', 'DRB1*03:09',
                 'DRB1*04:10', 'DRB5*01:01']
    __supported_length = [9]
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

    __alleles = ['B*27:20', 'B*83:01', 'A*32:15', 'B*15:17', 'B*40:13', 'A*24:02', 'A*24:03', 'B*53:01', 'B*15:01',
                 'B*27:05', 'B*42:01', 'B*39:01', 'B*38:01', 'A*23:01', 'A*25:01', 'C*04:01', 'A*29:02', 'A*02:06',
                 'A*02:01', 'A*02:02', 'A*02:03', 'A*26:02', 'A*26:03', 'A*26:01', 'C*03:03', 'E*01:01', 'E*01:03',
                 'B*58:01', 'A*31:01', 'C*06:02', 'B*07:02', 'A*66:01', 'B*57:01', 'A*68:01', 'A*68:02', 'C*14:02',
                 'B*35:01', 'B*15:09', 'B*35:03', 'A*80:01', 'B*15:03', 'B*15:02', 'A*32:01', 'A*02:50', 'A*32:07',
                 'B*58:02', 'A*69:01', 'A*68:23', 'A*11:01', 'A*03:01', 'B*73:01', 'B*40:01', 'B*44:03', 'B*46:01',
                 'B*40:02', 'C*12:03', 'B*44:02', 'A*30:01', 'A*02:19', 'A*30:02', 'A*02:17', 'A*02:16', 'B*51:01',
                 'B*45:01', 'A*02:12', 'A*02:11', 'B*54:01', 'B*08:01', 'B*18:01', 'B*08:03', 'B*08:02', 'C*05:01',
                 'C*15:02', 'A*33:01', 'B*14:02', 'C*07:01', 'B*48:01', 'B*15:42', 'C*07:02', 'A*01:01', 'C*08:02']
    __supported_length = [8, 9, 10, 11]
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
        return super(SMM, self).predict(peptides, alleles=alleles, **kwargs).applymap(lambda x: math.pow(10, x))


class SMMPMBEC(APSSMEpitopePredictor):
    """
    Implements IEDBs SMMPMBEC PSSM method
    """

    __alleles = ['A*01:01', 'A*02:01', 'A*02:02', 'A*02:03', 'A*02:06', 'A*02:11', 'A*02:12', 'A*02:16', 'A*02:17',
                 'A*02:19', 'A*02:50', 'A*03:01', 'A*11:01', 'A*23:01', 'A*24:02', 'A*24:03', 'A*25:01', 'A*26:01',
                 'A*26:02', 'A*26:03', 'A*29:02', 'A*30:01', 'A*30:02', 'A*31:01', 'A*32:01', 'A*32:07', 'A*32:15',
                 'A*33:01', 'A*66:01', 'A*68:01', 'A*68:02', 'A*68:23', 'A*69:01', 'A*80:01', 'B*07:02', 'B*08:01',
                 'B*08:02', 'B*08:03', 'B*14:02', 'B*15:01', 'B*15:02', 'B*15:03', 'B*15:09', 'B*15:17', 'B*15:42',
                 'B*18:01', 'B*27:03', 'B*27:05', 'B*27:20', 'B*35:01', 'B*35:03', 'B*38:01', 'B*39:01', 'B*40:01',
                 'B*40:02', 'B*40:13', 'B*42:01', 'B*44:02', 'B*44:03', 'B*45:01', 'B*45:06', 'B*46:01', 'B*48:01',
                 'B*51:01', 'B*53:01', 'B*54:01', 'B*57:01', 'B*58:01', 'B*58:02', 'B*73:01', 'B*83:01', 'C*03:03',
                 'C*04:01', 'C*05:01', 'C*06:02', 'C*07:01', 'C*07:02', 'C*08:02', 'C*12:03', 'C*14:02', 'C*15:02',
                 'E*01:01', 'E*01:03']
    __supported_length = [8, 9, 10, 11]
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
        return super(SMMPMBEC, self).predict(peptides, alleles=alleles, **kwargs).applymap(lambda x: math.pow(10, x))


class ARB(APSSMEpitopePredictor):
    """
    Implements IEDBs ARB method
    """

    __alleles = ['A*01:01', 'A*02:01', 'A*02:02', 'A*02:03', 'A*02:06', 'A*02:11', 'A*02:12', 'A*02:16', 'A*02:19',
                 'A*02:50', 'A*03:01', 'A*11:01', 'A*23:01', 'A*24:02', 'A*24:03', 'A*25:01', 'A*26:01', 'A*26:02',
                 'A*26:03', 'A*29:02', 'A*30:01', 'A*30:02', 'A*31:01', 'A*32:01', 'A*33:01', 'A*68:01', 'A*68:02',
                 'A*69:01', 'A*80:01', 'B*07:02', 'B*08:01', 'B*08:02', 'B*08:03', 'B*15:01', 'B*15:02', 'B*15:03',
                 'B*15:09', 'B*15:17', 'B*18:01', 'B*27:03', 'B*27:05', 'B*35:01', 'B*38:01', 'B*39:01', 'B*40:01',
                 'B*40:02', 'B*44:02', 'B*44:03', 'B*45:01', 'B*46:01', 'B*48:01', 'B*51:01', 'B*53:01', 'B*54:01',
                 'B*57:01', 'B*58:01', 'B*73:01']
    __supported_length = [8, 9, 10, 11]
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