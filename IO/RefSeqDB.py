# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'walzer'

import MySQLdb
import csv
import re
import urllib2
import warnings
import logging

class RefSeqDB():
    def __init__(self, usr=None, host=None, pwd=None, db=None):
        """
        used to fetch sequences from given RefSeq id's either from BioMart if no credentials given else from a MySQLdb
        :param usr: db user e.g. = 'ucsc_annot_query'
        :param host: db host e.g. = "pride"
        :param pwd: pw for user e.g. = 'an0q3ry'
        :param db: db on host e.g. = "hg18_ucsc_annotation"
        """
        if usr and host and pwd and db:
            self.connection = MySQLdb.connect(user=usr, host=host, passwd=pwd, db=db)
        else:
            self.connection = None

        self.biomart_url = """http://biomart.org/biomart/martservice?query="""
        self.biomart_head = """
        <?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE Query>
            <Query client="true" processor="TSV" limit="-1" header="1" uniqueRows = "1" >
                <Dataset name="hsapiens_gene_ensembl" config="gene_ensembl_config">
        """.strip()
        self.biomart_tail = """
                </Dataset>
            </Query>
        """.strip()
        self.biomart_filter = """<Filter name="%s" value="%s" filter_list=""/>"""
        self.biomart_attribute = """<Attribute name="%s"/>"""

    def get_product_sequence(self, product_refseq):
        """
        fetches product sequence for the given id
        :param product_refseq: given refseq id
        :return: list of dictionaries of the requested sequence, the respective strand and the associated gene name
        """
        if self.connection:
            cur = self.connection.cursor()
            #query = "SELECT * FROM sbs_ncbi_mrna WHERE id='%s';"%('5')
            query = "SELECT refseq,seq FROM hg19_ucsc_annotation.refLink INNER JOIN sbs_ncbi_protein on protAcc=refseq WHERE mrnaAcc='%s';" % (
                product_refseq)  # gives NP plus prot seq
            #mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'SELECT refseq,seq FROM hg19_ucsc_annotation.refLink INNER JOIN sbs_ncbi_protein on protAcc=refseq WHERE mrnaAcc='NM_002486';'
            try:
                cur.execute(query)
                result = cur.fetchall()
                if result:
                    product_refseq = result[0][0]
                    product_sequence = result[0][1]
                else:
                    warnings.warn("An Error occured while executing query: " + query + "\n")
                    return None
            except MySQLdb.Error, message:
                warnings.warn(
                    "An Error occured while executing query: " + query + "\n" + "Error message: " + message[1])
                return None
            return [{product_refseq: product_sequence}]
        else:
            filter = None
            if product_refseq.startswith('NP_'):
                filter = "refseq_predicted"
            elif product_refseq.startswith('XP_'):
                filter = "refseq_peptide_predicted"
            else:
                warnings.warn("No correct transcript id: " + product_refseq)
                return None
            rq_n = self.biomart_head \
                + self.biomart_filter%(filter, str(product_refseq))  \
                + self.biomart_attribute%("peptide")  \
                + self.biomart_attribute%(filter)  \
                + self.biomart_attribute%("external_gene_id")  \
                + self.biomart_attribute%("strand")  \
                + self.biomart_tail

            tsvreader = csv.DictReader(urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
            return [x for x in tsvreader]
        return None

    def get_transcript_sequence(self, transcript_refseq):
        """
        fetches transcript sequence for the given id
        :param transcript_refseq:
        :return: list of dictionary of the requested sequence, the respective strand and the associated gene name
        """
        if self.connection:
            cur = self.connection.cursor()
            #query = "SELECT * FROM sbs_ncbi_mrna WHERE id='%s';"%('5')
            query = "SELECT refseq,seq FROM hg19_ucsc_annotation.sbs_ncbi_mrna WHERE refseq LIKE '" + transcript_refseq + "%'"
            try:
                cur.execute(query)
                result = cur.fetchall()
                if result:
                    transcript_refseq = result[0][0]
                    transcript_sequence = result[0][1]
                    if len(result) > 1:
                        warnings.warn("Ambiguous transcript refseq query: " + transcript_refseq + "\n")
                else:
                    warnings.warn("An Error occured while executing query: " + query + "\n")
                    return None
            except MySQLdb.Error, message:
                warnings.warn("An Error occured while executing query: " + query + "\n" + "Error message: " + message[1])
                return None
            return (transcript_refseq, transcript_sequence)
        else:
            filter = None
            if transcript_refseq.startswith('NM_'):
                filter = "refseq_mrna"
            elif transcript_refseq.startswith('XM_'):
                filter = "refseq_mrna_predicted"
            else:
                warnings.warn("No correct transcript id: " + transcript_refseq)
                return None
            rq_n = self.biomart_head \
                + self.biomart_filter%(filter, str(transcript_refseq))  \
                + self.biomart_attribute%("peptide")  \
                + self.biomart_attribute%(filter)  \
                + self.biomart_attribute%("external_gene_id")  \
                + self.biomart_attribute%("strand")  \
                + self.biomart_tail

            tsvreader = csv.DictReader(urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
            return [x for x in tsvreader]
        return None

    def get_variant_gene(self, chrom, start, stop):
        """
        fetches the important db ids and names for given chromosomal location
        :param chrom: integer value of the chromosome in question
        :param start: integer value of the variation start position on given chromosome
        :param stop: integer value of the variation stop position on given chromosome
        :return: The respective gene name, i.e. the first one reported
        """
        rq_n = self.biomart_head \
            + self.biomart_filter%("chromosome_name", str(chrom))  \
            + self.biomart_filter%("start", str(start))  \
            + self.biomart_filter%("end", str(stop))  \
            + self.biomart_attribute%("uniprot_genename")  \
            + self.biomart_tail

        tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        return tsvselect[0]['UniProt Gene Name']

    def get_variant_ids(self, chrom, start, stop):
        """
        fetches the important db ids and names for given chromosomal location
        :param chrom: integer value of the chromosome in question
        :param start: integer value of the variation start position on given chromosome
        :param stop: integer value of the variation stop position on given chromosome
        :return: The list of dicts of entries with transcript and protein ids (either NM+NP or XM+XP)
        """
        # TODO add alternative Filter name = "uniprot_genename" and implement 'proxy'
        rq_n = self.biomart_head \
            + self.biomart_filter%("chromosome_name", str(chrom))  \
            + self.biomart_filter%("start", str(start))  \
            + self.biomart_filter%("end", str(stop))  \
            + self.biomart_attribute%("ensembl_gene_id")  \
            + self.biomart_attribute%("ensembl_peptide_id")  \
            + self.biomart_attribute%("ensembl_transcript_id")  \
            + self.biomart_attribute%("strand")  \
            + self.biomart_attribute%("refseq_mrna")  \
            + self.biomart_attribute%("refseq_peptide")  \
            + self.biomart_attribute%("uniprot_swissprot")  \
            + self.biomart_tail

        tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
        result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                  if x['RefSeq Protein ID [e.g. NP_001005353]'] and x['RefSeq mRNA [e.g. NM_001195597]']}

        rq_x = self.biomart_head \
            + self.biomart_filter%("chromosome_name", str(chrom))  \
            + self.biomart_filter%("start", str(start))  \
            + self.biomart_filter%("end", str(stop))  \
            + self.biomart_attribute%("ensembl_gene_id")  \
            + self.biomart_attribute%("ensembl_peptide_id")  \
            + self.biomart_attribute%("ensembl_transcript_id")  \
            + self.biomart_attribute%("refseq_peptide_predicted")  \
            + self.biomart_attribute%("refseq_mrna_predicted")  \
            + self.biomart_tail

        tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_x)).read()).splitlines(), dialect='excel-tab')

        result2 = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                   if x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and x['RefSeq mRNA predicted [e.g. XM_001125684]']}

        result.update(result2)

        return result.values()
