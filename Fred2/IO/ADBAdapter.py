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

Transcript sequences are always only coding region! First coding base has index i=1
Variant position start at conding region of the transcript
Strands are always encoded with + and -
'''

EAdapterFields = (lambda **enums: type('Enum', (), enums))(GENE=0, STRAND=1, SEQ=2, TRANSID=3, PROTID=4)


class ADBAdapter:

    __metaclass__ = ABCMeta

    @abstractmethod
    def get_product_sequence(self, **kwargs):
        """
        fetches the product sequence for the given id
        :keyword unknown: given product id, implement other keywords as well if you DBAdapter is supposed to
        support different id types
        :return: the requested sequence
        """
        pass

    @abstractmethod
    def get_transcript_sequence(self, **kwargs):
        """
        Fetches transcript sequence for the given id
        :keyword unknown: given transcript id, implement other keywords as well if you DBAdapter is supposed to
        support different id types
        :return: the requested sequence
        """
        pass

    @abstractmethod
    def get_transcript_information(self, **kwargs):
        """
        Fetches transcript sequence for the given id
        :keyword unknown: given transcript id, implement other keywords as well if you DBAdapter is supposed to
        support different id types
        :return: list of dictionary of the requested sequence, the respective strand and the associated gene name
        """
        pass