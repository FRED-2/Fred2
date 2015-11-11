# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: IO.ADBAdapter
   :synopsis: Base class for specific DB-Adapters
.. moduleauthor:: schubert
"""
from abc import ABCMeta, abstractmethod


'''
:Contract:

Transcript sequences are always only coding region! First coding base has index i=0
Variant position start at conding region of the transcript
Strands are always encoded with + and -
'''

EAdapterFields = (lambda **enums: type('Enum', (), enums))(GENE=0, STRAND=1, SEQ=2, TRANSID=3, PROTID=4)
EIdentifierTypes = (lambda **enums: type('Enum', (), enums))(ENSEMBL=0, REFSEQ=1, PREDREFSEQ=2, UNIPROT=3, GENENAME=4, HGNC=4)


class ADBAdapter:

    __metaclass__ = ABCMeta

    @abstractmethod
    def get_product_sequence(self, product_id, **kwargs):
        """
        fetches the product sequence for the given id
        :keyword type: given id, is in the form of this type,found in EIdentifierTypes. It is to be documented if an ADBAdapter implementation overrides these types.
        :return: the requested sequence

        """
        pass

    @abstractmethod
    def get_transcript_sequence(self, transcript_id, **kwargs):
        """
        Fetches transcript sequence for the given id

        :keyword type: given id, is in the form of this type,found in EIdentifierTypes. It is to be documented if an ADBAdapter implementation overrides these types.
        :return: the requested sequence

        """
        pass

    @abstractmethod
    def get_transcript_information(self, transcript_id, **kwargs):
        """
        Fetches transcript sequence for the given id
        :keyword type: given id, is in the form of this type,found in EIdentifierTypes. It is to be documented if an ADBAdapter implementation overrides these types.
        :return: list of dictionary of the requested sequence, the respective strand and the associated gene name
        """
        pass