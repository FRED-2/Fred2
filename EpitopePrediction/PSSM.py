# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'walzer', 'schubert'

import collections, itertools, warnings,pandas
from Fred2.Core.Result import Result
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
            return getattr( __import__("Fred2.Data.PSSMMatrices", fromlist=[allele_model]), allele_model)


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

            #dynamicaly import prediction PSSMS for alleles and predict
            if length not in self.supportedLength:
                warnings.warn("Peptide length of %i not supported"%length, RuntimeWarning)
                continue

            for a in allales_string.iterkeys():
                try:
                    pssm = __load_allele_model(a, length)
                except Exception:
                    warnings.warn("No model found for %s with length %i"%(allales_string[a], length))
                    continue
                result[allales_string[a]] = {}
                ##here is the prediction and result object missing##
                for p in peps:
                    score = sum(pssm[i][p[i]] for i in xrange(length))
                    result[allales_string[a]][str(p)] = score

        df_result = Result.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i,self.name)) for i in df_result.index],
                                                        names=['Seq','Method'])
        return df_result


class Syfpeithi(APSSMEpitopePredictor):
    """
        Represents the Syfpeithi PSSM predictor
    """
    __alleles = ["A*03:01"]
    __supported_length = [9, 10, 11]
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
        pass

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