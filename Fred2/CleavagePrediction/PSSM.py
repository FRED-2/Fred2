# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'schubert'

import collections
import itertools
import pandas
import warnings
import abc

from Bio.Alphabet import IUPAC

from Fred2.Core.Base import ACleavagePrediction
from Fred2.Core.Protein import Protein
from Fred2.Core.Result import CleavagePredictionResult

'''
NOTE: This implementation only supports cleavage site prediction not fragment prediction
'''


class APSSMCleavagePredictor(ACleavagePrediction):
    """
        Abstract base class for PSSM predictions.

        Implements predict functionality

    """

    @abc.abstractproperty
    def cleavagePos(self):
        """
        parameter specifying the position of aa (within the prediction window) after which the sequence is cleaved

        :return:
        """
        raise NotImplementedError

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
            model = "%s_%s_%i"%(self.name, length)

            #TODO: what if there exists no allele model for this length?
            return getattr( __import__("Fred2.Data.CleaveagePSSMMatrices", fromlist=[model]), model)

        if isinstance(aa_seq, collections.Iterable):
            if any(not isinstance(p.alphabet, IUPAC.IUPACProtein) for p in aa_seq):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(p):p for p in aa_seq}
        else:
            if not isinstance(aa_seq.alphabet, IUPAC.IUPACProtein):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(aa_seq):aa_seq}


        length = min(self.supportedLength) if length is None else length
        if length not in self.supportedLength:
            raise ValueError("Length %i is not supported by this method"%length)

        #group peptides by length and
        result = {}
        result["Cleavage Score"] = {}
        try:
            pssm = __load_model(length)
        except ImportError:
            raise KeyError("No model found for %s with length %i"%(self.name, length))

        diff = length - self.cleavagePos
        for seq in pep_seqs.iterkeys():

            #convention for peptides its always the first transcript ID after sorting
            seq_id = pep_seqs[seq].transcript_id if isinstance(pep_seqs[seq], Protein) else \
                                                    pep_seqs[seq].transcript_ids.keys().sort()[0]

            #dynamicaly import prediction PSSMS for alleles and predict
            if len(seq) < length:
                warnings.warn("Sequence length of %i is to small for specified window of %i"%(len(seq),length), RuntimeWarning)
                continue

            for i in xrange(len(seq)+diff):
                if i < length:
                    score = 0.0
                elif i > (len(seq) - (length-1)):
                    score = 0.0
                else:
                    score = sum(pssm[i][aa] for i,aa in enumerate(seq[i:i+length]))

                result[seq[i - diff]][seq_id] = score
                    #print a, score, result

        df_result = CleavagePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i,self.name)) for i in df_result.index],
                                                        names=['ID','Method'])
        return df_result


class PCM(APSSMCleavagePredictor):
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
    def name(self):
        return self.__name

    def predict(self, peptides, **kwargs):
        super(PCM, self).predict(peptides, **kwargs)