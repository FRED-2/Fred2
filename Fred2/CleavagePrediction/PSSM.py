# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: CleavagePrediction.PSSM
   :synopsis: PSSM-based cleavage prediction methods.
.. moduleauthor:: schubert

"""

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
        This implementation only supports cleavage site prediction not fragment prediction.
        Implements predict functionality.
    """

    def predict(self, aa_seq, length=None, **kwargs):
        """
        Returns predictions for given peptides.

        :param aa_seq: A single :class:`~Fred2.Core.Peptide.Peptide` or `~Fred2.Core.Protein.Protein` or a list of
                       :class:`~Fred2.Core.Peptide` or :class:`~Fred2.Core.Protein.Protein`
        :type aa_seq: list(:class:`~Fred2.Core.Peptide.Peptide` or :class:`~Fred2.Core.Protein.Protein`)
                      or :class:`~Fred2.Core.Peptide`/:class:`~Fred2.Core.Protein.Protein`
        :param int length: The peptide length of the cleavage site model. If None the default value is used.
        :return: Returns a :class:`~Fred2.Core.Result.CleavageSitePredictionResult` object
        :rtype: :class:`~Fred2.Core.Result.CleavageSitePredictionResult`
        """
        def __load_model(length):
            model = "%s_%i"%(self.name, length)
            return getattr(__import__("Fred2.Data.pssms."+self.name+".mat."+model, fromlist=[model]),
                           model)

        if isinstance(aa_seq, Peptide) or isinstance(aa_seq, Protein):
            pep_seqs = {str(aa_seq): aa_seq}
        else:
            pep_seqs = {}
            for p in aa_seq:
                if not isinstance(p, Peptide) and not isinstance(p, Protein):
                    raise ValueError("Input is not of type Protein or Peptide")
                pep_seqs[str(p)] = p

        length = min(self.supportedLength) if length is None else length
        if length not in self.supportedLength:
            raise ValueError("Length %i is not supported by %s"%(length, self.name))

        #group peptides by length and
        result = {"Seq": {}, self.name: {}}

        try:
            pssm = __load_model(length)
        except ImportError:
            raise KeyError("No model found for %s with length %i"%(self.name, length))

        diff = length - self.cleavagePos
        for j,seq in enumerate(pep_seqs.keys()):

            seq_id = "seq_%i"%j
            p = pep_seqs[seq]

            if isinstance(p, Protein):
                if p.transcript_id:
                    seq_id = p.transcript_id
            else:
                for t in p.proteins.keys():
                    if t:
                        seq_id = t
                        break

            #dynamicaly import prediction PSSMS for alleles and predict
            if len(seq) < length:
                warnings.warn("Sequence length of %i is to small for specified window of %i"%(len(seq),length), RuntimeWarning)
                continue

            for i in range(len(seq)):
                if i < (length-1):

                    result["Seq"][(seq_id, i)] = seq[i]
                    result[self.name][(seq_id, i)] = 0.0
                else:
                    result[self.name][(seq_id, i)] = 0.0
                    result["Seq"][(seq_id, i)] = seq[i]
                    score = sum(
                        pssm.get(i, {}).get(aa, 0) for i, aa in enumerate(seq[i - (length - 1):(i + 1)])) + pssm.get(-1,
                        {}).get("con", 0)
                    result[self.name][(seq_id, i-diff)] = score

        if not result["Seq"]:
            raise ValueError("No predictions could be made for the given input.")
        df_result = CleavageSitePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i, j)) for i, j in df_result.index],
                                                        names=['ID', 'Pos'])
        return df_result


class PCM(APSSMCleavageSitePredictor):
    """
    Implements the PCM cleavage prediction method.

    .. note::

        Doennes, P., and Kohlbacher, O. (2005). Integrated modeling of the major events in the MHC class
        I antigen processing pathway. Protein Science, 14(8), 2132-2140.
    """

    __supported_length = frozenset([6])
    __name = "pcm"
    __cleavage_pos = 4
    __version = "1.0"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    @property
    def supportedLength(self):
        """A list of supported peptide lengths"""
        return self.__supported_length

    @property
    def cleavagePos(self):
        """
        Parameter specifying the position of aa (within the prediction window) after which the sequence is cleaved
        """
        return self.__cleavage_pos

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    def predict(self, peptides, length=None, **kwargs):
        """
        Returns predictions for given peptides.

        :param aa_seq: A single :class:`~Fred2.Core.Peptide.Peptide`/`~Fred2.Core.Protein.Protein` or a list of
                       :class:`~Fred2.Core.Peptide`/:class:`~Fred2.Core.Protein.Protein`
        :type aa_seq: list(:class:`~Fred2.Core.Peptide.Peptide` or :class:`~Fred2.Core.Protein.Protein`)
                      or :class:`~Fred2.Core.Peptide`/:class:`~Fred2.Core.Protein.Protein`
        :param int length: The peptide length of the cleavage site model. If None the default value is used.
        :return: Returns a :class:`~Fred2.Core.Result.CleavageSitePredictionResult` object
        :rtype: :class:`~Fred2.Core.Result.CleavageSitePredictionResult`
        """
        return super(PCM, self).predict(peptides, lenght=length, **kwargs)


class ProteaSMMConsecutive(APSSMCleavageSitePredictor):
    """
    Implements the ProteaSMM cleavage prediction method.

    .. note::
        Tenzer, S., et al. "Modeling the MHC class I pathway by combining predictions of proteasomal cleavage,
        TAP transport and MHC class I binding." Cellular and Molecular Life Sciences CMLS 62.9 (2005): 1025-1037.

    This model represents the consecutive proteasom

    The matrices are generated not using the preon-dataset since a recent study has show that including those
    worsened the results.

    .. note::

        Calis, Jorg JA, et al. "Role of peptide processing predictions in T cell epitope identification:
        contribution of different prediction programs." Immunogenetics (2014): 1-9.

    """

    __supported_length = frozenset([10])
    __name = "proteasmm_c"
    __cleavage_pos = 6
    __version = "1.0"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    @property
    def supportedLength(self):
        """A list of supported peptide lengths"""
        return self.__supported_length

    @property
    def cleavagePos(self):
        """
        Parameter specifying the position of aa (within the prediction window) after which the sequence is cleaved
        """
        return self.__cleavage_pos

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    def predict(self, peptides, length=None, **kwargs):
        """
        Returns predictions for given peptides.

        :param aa_seq: A single :class:`~Fred2.Core.Peptide.Peptide`/`~Fred2.Core.Protein.Protein` or a list of
                       :class:`~Fred2.Core.Peptide`/:class:`~Fred2.Core.Protein.Protein`
        :type aa_seq: list(:class:`~Fred2.Core.Peptide.Peptide` or :class:`~Fred2.Core.Protein.Protein`)
                      or :class:`~Fred2.Core.Peptide`/:class:`~Fred2.Core.Protein.Protein`
        :param int length: The peptide length of the cleavage site model. If None the default value is used.
        :return: Returns a :class:`~Fred2.Core.Result.CleavageSitePredictionResult` object
        :rtype: :class:`~Fred2.Core.Result.CleavageSitePredictionResult`
        """
        return super(ProteaSMMConsecutive, self).predict(peptides, length=length, **kwargs)


class ProteaSMMImmuno(APSSMCleavageSitePredictor):
    """
    Implements the ProteaSMM cleavage prediction method.

    .. note::

        Tenzer, S., et al. "Modeling the MHC class I pathway by combining predictions of proteasomal cleavage,
        TAP transport and MHC class I binding." Cellular and Molecular Life Sciences CMLS 62.9 (2005): 1025-1037.

    This model represents the immuno proteasom

    The matrices are generated not using the preon-dataset since a recent study has show that including those
    worsened the results.

    .. note::

        Calis, Jorg JA, et al. "Role of peptide processing predictions in T cell epitope identification:
        contribution of different prediction programs." Immunogenetics (2014): 1-9.

    """

    __supported_length = frozenset([10])
    __name = "proteasmm_i"
    __cleavage_pos = 6
    __version = "1.0"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    @property
    def supportedLength(self):
        """A list of supported peptide lengths"""
        return self.__supported_length

    @property
    def cleavagePos(self):
        """
        Parameter specifying the position of aa (within the prediction window) after which the sequence is cleaved
        """
        return self.__cleavage_pos

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    def predict(self, peptides, length=None, **kwargs):
        """
        Returns predictions for given peptides.

        :param aa_seq: A single :class:`~Fred2.Core.Peptide.Peptide`/`~Fred2.Core.Protein.Protein` or a list of
                       :class:`~Fred2.Core.Peptide`/:class:`~Fred2.Core.Protein.Protein`
        :type aa_seq: list(:class:`~Fred2.Core.Peptide.Peptide` or :class:`~Fred2.Core.Protein.Protein`)
                      or :class:`~Fred2.Core.Peptide`/:class:`~Fred2.Core.Protein.Protein`
        :param int length: The peptide length of the cleavage site model. If None the default value is used.
        :return: Returns a :class:`~Fred2.Core.Result.CleavageSitePredictionResult` object
        :rtype: :class:`~Fred2.Core.Result.CleavageSitePredictionResult`
        """
        return super(ProteaSMMImmuno, self).predict(peptides, length=length, **kwargs)


class APSSMCleavageFragmentPredictor(ACleavageFragmentPrediction):
    """
        Abstract base class for PSSM predictions.

        This implementation only supports cleavage fragment prediction not site prediction

        Implements predict functionality
    """

    @abc.abstractproperty
    def trailingN(self):
        """
        The number of trailing residues at the N-terminal of the peptide used for prediction
        """
        raise NotImplementedError

    @abc.abstractproperty
    def tralingC(self):
        """
        The number of trailing residues at the C-terminal of the peptide used for prediction
        """
        raise NotImplementedError

    def predict(self, peptides,  **kwargs):
        """
        Takes peptides plus their trailing C and N-terminal residues to predict
        the probability that this n-mer was produced by proteasomal cleavage. It returns the score and
        the peptide sequence in a AResult object. Row-IDs are the peitopes column is the prediction score.

        :param peptides: A list of peptide objects or a single peptide object
        :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`) or :class:`~Fred2.Core.Peptide.Peptide`
        :return: Returns a :class:`Fred2.Core.Result.CleavageFragmentPredictionResult` object
        :rtype: :class:`Fred2.Core.Result.CleavageFragmentPredictionResult`
        """
        def __load_model(length):
            allele_model = "%%s_%i"%(self.name, length)
            return getattr(__import__("Fred2.Data.pssms."+self.name+".mat."+allele_model, fromlist=[allele_model]),
                           allele_model)

        if isinstance(peptides, Peptide):
            pep_seqs = {str(peptides):peptides}
        else:
            pep_seqs = {}
            for p in peptides:
                if not isinstance(p, Peptide):
                    raise ValueError("Input is not of type Protein or Peptide")
                pep_seqs[str(p)] = p

        result = {self.name:{}}
        for length, peps in itertools.groupby(iter(pep_seqs.keys()), key= lambda x: len(x)):
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
                score = sum(pssm[i][aa] for i, aa in enumerate(p))+pssm.get(-1, {}).get("con", 0)
                pep = pep_seqs[p][self.trailingN:-self.tralingC]
                result[self.name][pep] = score

        if not result:
            raise ValueError("No predictions could be made for the given input.")
        df_result = CleavageFragmentPredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i, self.name)) for i in df_result.index],
                                                        names=['Seq', 'Method'])
        return df_result


class PSSMGinodi(APSSMCleavageFragmentPredictor):
    """
    Implements the Cleavage Fragment prediction method of Ginodi et al.

    .. note::

        Ido Ginodi, Tal Vider-Shalit, Lea Tsaban, and Yoram Louzoun
        Precise score for the prediction of peptides cleaved by the proteasome
        Bioinformatics (2008) 24 (4): 477-483

    """

    __supported_lengths=[11]
    __name="ginodi"
    __trailingN = 1
    __trailingC = 1
    __version = "1.0"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    @property
    def cleavagePos(self):
        """
        Parameter specifying the position of aa (within the prediction window) after which the sequence is cleaved
        """
        return 0

    @property
    def supportedLength(self):
        """A list of supported peptide lengths"""
        return self.__supported_lengths

    @property
    def trailingN(self):
        """
        The number of trailing residues at the N-terminal of the peptide used for prediction
        """
        return self.__trailingN

    @property
    def tralingC(self):
        """
        The number of trailing residues at the C-terminal of the peptide used for prediction
        """
        return self.__trailingC

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    def predict(self, peptides,  **kwargs):
        """
        Takes peptides plus their trailing C and N-terminal residues to predict
        the probability that this n-mer was produced by proteasomal cleavage. It returns the score and
        the peptide sequence in a AResult object. Row-IDs are the peitopes column is the prediction score.

        :param peptides: A list of peptide objects or a single peptide object
        :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`) or :class:`~Fred2.Core.Peptide.Peptide`
        :return: Returns a :class:`Fred2.Core.Result.CleavageFragmentPredictionResult` object
        :rtype: :class:`Fred2.Core.Result.CleavageFragmentPredictionResult`
        """
        def __load_model(length):
            allele_model = "%s_%i"%(self.name, length)
            return getattr(__import__("Fred2.Data.pssms."+self.name+".mat."+allele_model, fromlist=[allele_model]),
                           allele_model)

        if isinstance(peptides, Peptide):
            pep_seqs = {str(peptides): peptides}
        else:
            pep_seqs = {}
            for p in peptides:
                if not isinstance(p, Peptide) and not isinstance(p, Protein):
                    raise ValueError("Input is not of type Protein or Peptide")
                pep_seqs[str(p)] = p

        result = {self.name: {}}
        for length, peps in itertools.groupby(iter(pep_seqs.keys()), key=lambda x: len(x)):
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
                score = pssm[0][p[0]]+pssm[1][p[1]] + sum(pssm[2][aa] for aa in p[2:-2])\
                        + pssm[3][p[-2]] + pssm[4][p[-1]]
                pep = pep_seqs[p][self.trailingN: -self.tralingC]
                result[self.name][pep] = score

        if not result:
            raise ValueError("No predictions could be made for the given input.")
        df_result = CleavageFragmentPredictionResult.from_dict(result)
        return df_result