# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
from abc import ABCMeta, abstractmethod

__author__ = 'schubert'

'''
:Contract:

Transcript sequences are always only coding region! First coding base has index i=1
Variant position start at conding region of the transcript
Strands are always encoded with + and -
'''

EAdapterFields = (lambda **enums: type('Enum', (), enums))(GENE=0, STRAND=1, SEQ=2, TRANSID=3, PROTID=4)



class ADBAdapter:
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_product_sequence(self, product_refseq):
        """
        fetches product sequence for the given id

        :param product_refseq: given refseq id
        :return: list of dictionaries of the requested sequence, the respective strand and the associated gene name
        """
        pass

    @abstractmethod
    def get_transcript_sequence(self, transcript_refseq):
        """
        Fetches transcript sequence for the given id
        
        :param transcript_refseq:
        :return: list of dictionary of the requested sequence, the respective strand and the associated gene name
        """
        pass

    @abstractmethod
    def get_transcript_information(self, transcript_refseq):
        pass