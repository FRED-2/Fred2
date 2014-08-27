from abc import ABCMeta, abstractmethod

__author__ = 'schubert'


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
