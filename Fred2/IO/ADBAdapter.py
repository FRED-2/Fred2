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


class ADBAdapter(metaclass=ABCMeta):

    @abstractmethod
    def get_product_sequence(self, product_id, **kwargs):
        """
        Fetches the product sequence for the given id

        :param str product_id: The product ID as string
        :keyword type: Given id, is in the form of this type,found in :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`. It
                       is to be documented if an ADBAdapter implementation overrides these types.
        :type type: :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`
        :return: The requested sequence
        :rtype: str
        """
        raise NotImplementedError

    @abstractmethod
    def get_transcript_sequence(self, transcript_id, **kwargs):
        """
        Fetches transcript sequence for the given id

        :param str transcript_id: The transcript ID as string
        :keyword type: Given id, is in the form of this type,found in :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`. It
                       is to be documented if an ADBAdapter implementation overrides these types.
        :type type: :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`
        :return: The requested sequence
        :rtype: str
        """
        raise NotImplementedError

    @abstractmethod
    def get_transcript_information(self, transcript_id, **kwargs):
        """
        Fetches transcript sequence for the given id

        :param str transcript_id: The transcript ID as string
        :keyword type: Given id, is in the form of this type,found in :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`. It
                       is to be documented if an ADBAdapter implementation overrides these types.
        :return: list of dictionary of the requested sequence, the respective strand and the associated gene name
        :rtype: list(dict)
        """
        raise NotImplementedError