# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: DummyAdaper
   :synopsis: Contains a pseudo data base adapter for testing purposes.
.. moduleauthor:: schubert, brachvogel
"""
import copy
from Fred2.IO.ADBAdapter import ADBAdapter, EAdapterFields


class DummyAdapter(ADBAdapter):

    def __init__(self):
        pass

    def get_product_sequence(self, product_refseq):
        # TODO: also implement this one?
        pass

    def get_transcript_sequence(self, transcript_refseq):
        # TODO: also implement this one?
        pass

    def get_transcript_information(self, transcript_refseq):
        """
        At the moment we only use this method.
        :param transcript_refseq: Refseq id of transcript
        :type transcript_refseq: str.
        :return: Dictionary with (EAdapterFields: <field content>
        relevant: GENE = gene id, STRAND = +/-, SEQ = transcript sequence
        """
        tsc_1 = {
            EAdapterFields.SEQ: "AAAAACCCCCGGGGG", # 15 * C
            EAdapterFields.GENE: "gene_1", # gene id
            EAdapterFields.STRAND: "+", # normal 5' to 3'
        }
        tsc_2 = {
            EAdapterFields.SEQ: "GGGGGCCCCCAAAAA", # 15 * C
            EAdapterFields.GENE: "gene_1", # gene id
            EAdapterFields.STRAND: "+", # normal 5' to 3'
        }

        res = {
            "tsc_1": tsc_1,
            "tsc_2": tsc_2
        }
        return copy.deepcopy(res[transcript_refseq])

