# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'schubert'

import collections
import itertools
import pandas
import warnings
import abc

from Fred2.Core.Base import ACleavageSitePrediction, ACleavageFragmentPrediction
from Fred2.Core.Protein import Protein
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Result import CleavageSitePredictionResult, CleavageFragmentPredictionResult


class APSSMCleavageSitePredictor(ACleavageSitePrediction):
    """
        Abstract base class for PSSM predictions.

        This implementation only supports cleavage site prediction not fragment prediction

        Implements predict functionality

    """

    def predict(self, aa_seq, length=None, **kwargs):
        """
        Returns predictions for given peptides an alleles. If no alleles are given, predictions for all available models
        are made.

        :param list(Peptide)/Peptide peptides: A single Peptide or a list of Peptides
        :param list(Alleles) alleles: a list of Alleles
        :param kwargs: optional parameter (not used yet)
        :return: Returns a Result object with the prediction results
        """
        def __load_model(length):
            model = "%s_%i"%(self.name, length)

            #TODO: what if there exists no allele model for this length?
            return getattr( __import__("Fred2.Data.CleaveagePSSMMatrices", fromlist=[model]), model)

        if isinstance(aa_seq, collections.Iterable):
            for p in aa_seq:
                print isinstance(p, Peptide)
            if any((not isinstance(p, Peptide)) and (not isinstance(p, Protein)) for p in aa_seq):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(p):p for p in aa_seq}
        else:
            if (not isinstance(aa_seq, Peptide)) or (not isinstance(aa_seq, Protein)):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(aa_seq):aa_seq}


        length = min(self.supportedLength) if length is None else length
        if length not in self.supportedLength:
            raise ValueError("Length %i is not supported by this method"%length)

        #group peptides by length and
        result = {"Seq":{},self.name:{}}

        try:
            pssm = __load_model(length)
        except ImportError:
            raise KeyError("No model found for %s with length %i"%(self.name, length))

        diff = length - self.cleavagePos
        for j,seq in enumerate(pep_seqs.iterkeys()):

            #convention for peptides its always the first transcript ID after sorting
            try:
                seq_id = pep_seqs[seq].transcript_id if isinstance(pep_seqs[seq], Protein) else \
                                                    pep_seqs[seq].transcripts.keys().sort()[0]
            except Exception:
                seq_id = "seq_%i"%j

            #dynamicaly import prediction PSSMS for alleles and predict
            if len(seq) < length:
                warnings.warn("Sequence length of %i is to small for specified window of %i"%(len(seq),length), RuntimeWarning)
                continue

            for i in xrange(len(seq)):
                if i < (length-1):

                    result["Seq"][(seq_id, i)] = seq[i]
                    result[self.name][(seq_id, i)] = 0.0
                else:
                    result[self.name][(seq_id, i)] = 0.0
                    result["Seq"][(seq_id, i)] = seq[i]
                    score = sum(pssm[i][aa] for i,aa in enumerate(seq[i-(length-1):(i+1)]))
                    result[self.name][(seq_id, i-diff)] = score

        df_result = CleavageSitePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i,j)) for i,j in df_result.index],
                                                        names=['ID','Pos'])

        return df_result


class PCM(APSSMCleavageSitePredictor):
    """
        Implements the PCM cleavage prediction method

        :input list(Peptide) peptides: list of peptides for which cleave probability prediction should be performed
    """

    __supported_length = [6]
    __name = "pcm"
    __cleavage_pos = 4

    @property
    def supportedLength(self):
        return self.__supported_length

    @property
    def cleavagePos(self):
        return self.__cleavage_pos

    @property
    def name(self):
        return self.__name

    def predict(self, peptides, **kwargs):
        return super(PCM, self).predict(peptides, **kwargs)


class APSSMCleavageFragmentPredictor(ACleavageFragmentPrediction):
    """
        Abstract base class for PSSM predictions.

        This implementation only supports cleavage fragment prediction not site prediction

        Implements predict functionality
    """

    @abc.abstractproperty
    def trailingN(self):
        """
        returns the number of trailing residues at the N terminal of the peptide used for prediction

        :return:
        """
        raise NotImplementedError

    @abc.abstractproperty
    def tralingC(self):
        """
        returns the number of trailing residues at the N terminal of the peptide used for prediction

        :return:
        """
        raise NotImplementedError

    def predict(self, peptides,  **kwargs):
        """
        Takes peptides plus their trailing C and N-terminal residues to predict
        the probability that this 9mer was produced by proteasomal cleavage. It returns the score and
        the peptide sequence in a Result object. Row-IDs are the peitopes column is the prediction score.

        :param _aa_seq:
        :param kwargs:
        :return:
        """
        def __load_model(length):
            allele_model = "%%s_%i"%(self.name, length)

            #TODO: what if there exists no allele model for this length?
            return getattr( __import__("Fred2.Data.CleaveagePSSMMatrices", fromlist=[allele_model]), allele_model)

        if isinstance(peptides, collections.Iterable):
            pep_seqs = {str(p):p for p in peptides}
        else:
            pep_seqs = {str(peptides):peptides}

        result = {self.name:{}}
        for length, peps in itertools.groupby(pep_seqs.iterkeys(), key= lambda x: len(x)):
            peps = list(peps)
            #dynamicaly import prediction PSSMS for alleles and predict
            if length not in self.supportedLength:
                warnings.warn("Peptide length of %i not supported"%length, RuntimeWarning)
                continue

            #load pssm matrices
            try:
                pssm = __load_model(length)
            except ImportError:
                raise KeyError("No model found for %s with length %i"%(self.name, length))

            for p in peps:
                score = sum(pssm[i][aa] for i,aa in enumerate(p))
                pep = pep_seqs[p][self.trailingN:-self.tralingC]
                result[self.name][pep] = score

        df_result = CleavageFragmentPredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i,self.name)) for i in df_result.index],
                                                        names=['Seq','Method'])
        return df_result


class PSSMGinodi(APSSMCleavageFragmentPredictor):
    """
        Implements the Cleavage Fragment prediction method of Ginodi et al.

        Ido Ginodi, Tal Vider-Shalit, Lea Tsaban, and Yoram Louzoun
        Precise score for the prediction of peptides cleaved by the proteasome
        Bioinformatics (2008) 24 (4): 477-483

    """

    __supported_lengths=[11]
    __name="ginodi"
    __trailingN = 1
    __trailingC = 1

    @property
    def cleavagePos(self):
        return 0

    @property
    def supportedLength(self):
        return self.__supported_lengths

    @property
    def trailingN(self):
        return self.__trailingN

    @property
    def tralingC(self):
        return self.__trailingC

    @property
    def name(self):
        return self.__name

    def predict(self, peptides,  **kwargs):
        def __load_model(allele):
            allele_model = "%s_%i"%(self.name, length)

            #TODO: what if there exists no allele model for this length?
            return getattr( __import__("Fred2.Data.CleaveagePSSMMatrices", fromlist=[allele_model]), allele_model)

        if isinstance(peptides, collections.Iterable):
            pep_seqs = {str(p):p for p in peptides}
        else:
            pep_seqs = {str(peptides):peptides}

        for length, peps in itertools.groupby(pep_seqs.iterkeys(), key= lambda x: len(x)):
            peps = list(peps)
            #dynamicaly import prediction PSSMS for alleles and predict
            if length not in self.supportedLength:
                warnings.warn("Peptide length of %i not supported"%length, RuntimeWarning)
                continue

            #load pssm matrices
            try:
                pssm = __load_model(length)
            except ImportError:
                raise KeyError("No model found for %s with length %i"%(self.name, length))

            result = {self.name:{}}
            for p in peps:
                score = pssm[0][p[0]]+pssm[1][p[1]] + sum(pssm[2][aa] for aa in p[2:-2])\
                        + pssm[3][p[-2]] + pssm[4][p[-1]]
                pep = pep_seqs[p][self.trailingN: -self.tralingC]
                result[self.name][pep] = score

        df_result = CleavageFragmentPredictionResult.from_dict(result)
        return df_result