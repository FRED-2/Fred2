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
        :return: tuple of id and transcript sequence as strings
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
            return (product_refseq, product_sequence)
        else:
            rq = "http://biomart.org/biomart/martservice?query=%3C?xml%20version=%221.0%22%20encoding=%22UTF-8%22?%3E%3C!DOCTYPE%20Query%3E%3CQuery%20virtualSchemaName%20=%20%22default%22%20formatter%20=%20%22FASTA%22%20header%20=%20%220%22%20uniqueRows%20=%20%220%22%20count%20=%20%22%22%20datasetConfigVersion%20=%20%220.6%22%20%3E%3CDataset%20name%20=%20%22hsapiens_gene_ensembl%22%20interface%20=%20%22default%22%20%3E%3CFilter%20name%20=%20%refseq_mrna%22%20value%20=%20%22@XXX@%22/%3E%3CAttribute%20name%20=%20%22peptide%22%20/%3E%3CAttribute%20name%20=%20%22ensembl_gene_id%22%20/%3E%20%3CAttribute%20name%20=%20%22ensembl_transcript_id%22%20/%3E%20%3CAttribute%20name%20=%20%22ensembl_peptide_id%22%20/%3E%20%3CAttribute%20name%20=%20%22refseq_mrna%22%20/%3E%20%3C/Dataset%3E%3C/Query%3E"
            urllib2.urlopen(rq).read()
            rq = rq.replace('@XXX@', product_refseq)
            return (product_refseq, ''.join(urllib2.urlopen(rq).read().split('\n')[1:]))
        return None

    def get_transcript_sequence(self, transcript_refseq):
        """
        fetches transcript sequence for the given id
        :param transcript_refseq:
        :return: tuple of id and product sequence as string
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
            rq = "http://biomart.org/biomart/martservice?query=%3C?xml%20version=%221.0%22%20encoding=%22UTF-8%22?%3E%3C!DOCTYPE%20Query%3E%3CQuery%20virtualSchemaName%20=%20%22default%22%20formatter%20=%20%22TSV%22%20header%20=%20%220%22%20uniqueRows%20=%20%220%22%20count%20=%20%22%22%20datasetConfigVersion%20=%20%220.6%22%20%3E%20%3CDataset%20name%20=%20%22hsapiens_gene_ensembl%22%20interface%20=%20%22default%22%20%3E%20%3CFilter%20name%20=%20%22refseq_mrna%22%20value%20=%20%22@XXX@%22/%3E%20%3CAttribute%20name%20=%20%22ensembl_gene_id%22%20/%3E%20%3CAttribute%20name%20=%20%22ensembl_transcript_id%22%20/%3E%20%3CAttribute%20name%20=%20%22ensembl_peptide_id%22%20/%3E%20%3CAttribute%20name%20=%20%22refseq_mrna%22%20/%3E%20%3C/Dataset%3E%3C/Query%3E"
            rq = rq.replace('@XXX@', transcript_refseq)
            return (transcript_refseq, ''.join(urllib2.urlopen(rq).read().split('\n')[1:]))
        return None

    def get_variant_gene(self, chrom, start, stop):
        """
        fetches the important db ids and names for given chromosomal location
        :param chrom: integer value of the chromosome in question
        :param start: integer value of the variation start position on given chromosome
        :param stop: integer value of the variation stop position on given chromosome
        :return: tuple of location and list of dicts
        """

        rq = "http://biomart.org/biomart/martservice?query=%3C?xml%20version=%221.0%22%20encoding=%22UTF-8%22?%3E%3C!DOCTYPE%20Query%3E%3CQuery%20virtualSchemaName%20=%20%22default%22%20formatter%20=%20%22TSV%22%20header%20=%20%220%22%20uniqueRows%20=%20%220%22%20count%20=%20%22%22%20datasetConfigVersion%20=%20%220.6%22%20%3E%3CDataset%20name%20=%20%22hsapiens_gene_ensembl%22%20interface%20=%20%22default%22%20%3E%3CFilter%20name%20=%20%22chromosome_name%22%20value%20=%20%22@chr@%22/%3E%3CFilter%20name%20=%20%22end%22%20value%20=%20%22@stop@%22/%3E%3CFilter%20name%20=%20%22start%22%20value%20=%20%22@start@%22/%3E%3CAttribute%20name%20=%20%22external_gene_id%22%20/%3E%3CAttribute%20name%20=%20%22ensembl_gene_id%22%20/%3E%3CAttribute%20name%20=%20%22ensembl_transcript_id%22%20/%3E%3CAttribute%20name%20=%20%22ensembl_peptide_id%22%20/%3E%3CAttribute%20name%20=%20%22refseq_mrna%22%20/%3E%3CAttribute%20name%20=%20%22refseq_peptide%22%20/%3E%3C/Dataset%3E%3C/Query%3E"
        rq = rq.replace('@chr@', str(chrom)).replace('@start@', str(start)).replace('@stop@', str(stop))
        header = "external_gene_id\tensembl_gene_id\tensembl_transcript_id\tensembl_peptide_id\trefseq_mrna\trefseq_peptide"
        tsvreader = csv.DictReader((header+'\n'+urllib2.urlopen(rq).read()).splitlines(), dialect='excel-tab')
        return ((chrom, start, stop), [x for x in tsvreader if x['refseq_mrna']])  # TODO make _1_ definitive response

    def get_variant_ids(self, chrom, start, stop):
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

        # print cou
        return result
