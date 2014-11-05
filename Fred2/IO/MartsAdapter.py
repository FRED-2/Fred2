# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'walzer', 'schubert'

import MySQLdb
import csv
import urllib2
import warnings
import logging

from Fred2.IO.ADBAdapter import ADBAdapter, EAdapterFields

class MartsAdapter(ADBAdapter):
    def __init__(self, usr=None, host=None, pwd=None, db=None):
        """
        used to fetch sequences from given RefSeq id's either from BioMart if no credentials given else from a MySQLdb
        :param usr: db user e.g. = 'ucsc_annot_query'
        :param host: db host e.g. = "pride"
        :param pwd: pw for user e.g. = 'an0q3ry'
        :param db: db on host e.g. = "hg18_ucsc_annotation"
        """
        self.ids_proxy = dict()
        self.gene_proxy = dict()
        self.sequence_proxy = dict()

        if usr and host and pwd and db:
            self.connection = MySQLdb.connect(user=usr, host=host, passwd=pwd, db=db)
        else:
            self.connection = None

        self.biomart_url = "http://biomart.org/biomart/martservice?query="
        self.new_biomart_url = """http://central.biomart.org/biomart/martservice?query="""
        self.biomart_head = """
        <?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE Query>
            <Query client="true" processor="TSV" limit="-1" header="1" uniqueRows = "1" >
                <Dataset name="%s" config="%s">
        """.strip()
        self.biomart_tail = """
                </Dataset>
            </Query>
        """.strip()
        self.biomart_filter = """<Filter name="%s" value="%s" filter_list=""/>"""
        self.biomart_attribute = """<Attribute name="%s"/>"""

    #TODO: refactor ... function based on old code
    def get_product_sequence(self, product_refseq, _db="hsapiens_gene_ensembl", _dataset='gene_ensembl_config'):
        """
        fetches product sequence for the given id
        :param product_refseq: given refseq id
        :return: list of dictionaries of the requested sequence, the respective strand and the associated gene name
        """

        if product_refseq in self.sequence_proxy:
            return self.sequence_proxy[product_refseq]

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
            self.sequence_proxy[product_refseq] = product_sequence
            return [{product_refseq: product_sequence}]
        else:
            filter = None
            if product_refseq.startswith('NP_'):
                filter = "refseq_peptide"
            elif product_refseq.startswith('XP_'):
                filter = "refseq_peptide_predicted"
            elif product_refseq.startswith('ENS'):
                filter = "ensembl_peptide_id"
            else:
                warnings.warn("No correct product id: " + product_refseq)
                return None
            rq_n = self.biomart_head%(_db,_dataset) \
                + self.biomart_filter%(filter, str(product_refseq))  \
                + self.biomart_attribute%("peptide")  \
                + self.biomart_attribute%(filter)  \
                + self.biomart_attribute%("external_gene_id")  \
                + self.biomart_attribute%("strand")  \
                + self.biomart_tail

            logging.warn(rq_n)

            tsvreader = csv.DictReader(urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
            tsvselect = [x for x in tsvreader]
            if not tsvselect:
                return None
            self.sequence_proxy[product_refseq] = tsvselect[0]["Protein"]
            return self.sequence_proxy[product_refseq]

    #TODO: refactor ... function based on old code
    def get_transcript_sequence(self, transcript_refseq, _db="hsapiens_gene_ensembl", _dataset='gene_ensembl_config'):
        """
        Fetches transcript sequence for the given id
        :param transcript_refseq:
        :return: list of dictionary of the requested sequence, the respective strand and the associated gene name
        """
        # TODO transcript_refseq to transcript_id, sniff which, query according
        if transcript_refseq in self.sequence_proxy:
            return self.sequence_proxy[transcript_refseq]

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
            self.sequence_proxy[transcript_refseq] = transcript_sequence
            return [{transcript_refseq: transcript_sequence}]
        else:
            filter = None
            if transcript_refseq.startswith('NM_'):
                filter = "refseq_mrna"
            elif transcript_refseq.startswith('XM_'):
                filter = "refseq_mrna_predicted"
            elif transcript_refseq.startswith('ENS'):
                filter = "ensembl_transcript_id"
            else:
                warnings.warn("No correct transcript id: " + transcript_refseq)
                return None
            rq_n = self.biomart_head%(_db, _dataset) \
                + self.biomart_filter%(filter, str(transcript_refseq))  \
                + self.biomart_attribute%(filter)  \
                + self.biomart_attribute%("coding")  \
                + self.biomart_attribute%("external_gene_id")  \
                + self.biomart_attribute%("strand")  \
                + self.biomart_tail

            tsvreader = csv.DictReader(urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
            tsvselect = [x for x in tsvreader]
            self.sequence_proxy[transcript_refseq] = tsvselect[0]['Coding sequence']
            return self.sequence_proxy[transcript_refseq]

    def get_transcript_information(self, transcript_refseq, _db="hsapiens_gene_ensembl", _dataset='gene_ensembl_config'):
        """
        It also already uses the Field-Enum for DBAdapters

        Fetches transcript sequence for the given id
        :param transcript_refseq:
        :return: list of dictionary of the requested sequence, the respective strand and the associated gene name
        """
        # TODO transcript_refseq to transcript_id, sniff which, query according
        if transcript_refseq in self.sequence_proxy:
            return self.sequence_proxy[transcript_refseq]

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
            self.ids_proxy[transcript_refseq] = transcript_sequence
            return [{transcript_refseq: transcript_sequence}]
        else:
            filter = None
            if transcript_refseq.startswith('NM_'):
                filter = "refseq_mrna"
            elif transcript_refseq.startswith('XM_'):
                filter = "refseq_mrna_predicted"
            elif transcript_refseq.startswith('ENS'):
                filter = "ensembl_transcript_id"
            else:
                warnings.warn("No correct transcript id: " + transcript_refseq)
                return None
            rq_n = self.biomart_head%(_db, _dataset) \
                + self.biomart_filter%(filter, str(transcript_refseq))  \
                + self.biomart_attribute%(filter)  \
                + self.biomart_attribute%("coding")  \
                + self.biomart_attribute%("strand")  \
                + self.biomart_attribute%("external_gene_id") \
                + self.biomart_tail

            print "Transcript information ",rq_n
            tsvreader = csv.DictReader(urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
            tsvselect = [x for x in tsvreader]
            if not tsvselect:
                warnings.warn("No entry found for ID %s"%transcript_refseq)
                return None

            self.ids_proxy[transcript_refseq] = {EAdapterFields.SEQ: tsvselect[0]['Coding sequence'],
                                                      EAdapterFields.GENE: tsvselect[0]['Associated Gene Name'],
                                                      EAdapterFields.STRAND: "-" if int(tsvselect[0]['Strand']) < 0
                                                      else "+"}
            return self.ids_proxy[transcript_refseq]

    #TODO: refactor ... function based on old code
    def get_variant_gene(self, chrom, start, stop, _db="hsapiens_gene_ensembl", _dataset='gene_ensembl_config'):
        """
        Fetches the important db ids and names for given chromosomal location
        :param chrom: integer value of the chromosome in question
        :param start: integer value of the variation start position on given chromosome
        :param stop: integer value of the variation stop position on given chromosome
        :return: The respective gene name, i.e. the first one reported
        """
        if chrom + start + stop in self.gene_proxy:
            return self.gene_proxy[chrom + start + stop]

        rq_n = self.biomart_head%(_db, _dataset) \
            + self.biomart_filter%("chromosome_name", str(chrom))  \
            + self.biomart_filter%("start", str(start))  \
            + self.biomart_filter%("end", str(stop))  \
            + self.biomart_attribute%("uniprot_genename")  \
            + self.biomart_tail

        tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        if tsvselect and tsvselect[0]:
            self.gene_proxy[chrom + start + stop] = tsvselect[0]['UniProt Gene Name']
            return tsvselect[0]['UniProt Gene Name']
        else:
            logging.warning(','.join([str(chrom), str(start), str(stop)]) + ' does not denote a known gene location')
            return ''

    def get_transcript_information_from_protein_id(self, **kwargs):
        """
        It also already uses the Field-Enum for DBAdapters

        Fetches transcript sequence for the given id
        :param transcript_refseq:
        :return: list of dictionary of the requested sequence, the respective strand and the associated gene name
        """
        _db = kwargs.get("_db","hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")
        filter = None
        db_id = ""
        if "refseq" in kwargs:
            filter = "refseq_peptide"
            db_id = kwargs["refseq"]
        elif "ensemble" in kwargs:
            filter = "ensembl_peptide_id"
            db_id = kwargs["ensemble"]
        elif "swiss_accid" in kwargs:
            filter = "uniprot_swissprot_accession"
            db_id = kwargs["swiss_accid"]
        elif "swiss_gene" in kwargs:
            filter= "uniprot_swissprot"
            db_id = kwargs["swiss_gene"]
        else:
            warnings.warn("No correct transcript id")
            return None
        rq_n = self.biomart_head%(_db, _dataset) \
               + self.biomart_filter%(filter, str(db_id)) \
               + self.biomart_attribute%(filter) \
               + self.biomart_attribute%("coding") \
               + self.biomart_attribute%("external_gene_id") \
               + self.biomart_attribute%("strand") \
               + self.biomart_tail

        tsvreader = csv.DictReader(urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        print tsvselect[0]
        if not tsvselect:
            warnings.warn("No entry found for ID %s"%db_id)
            return None

        self.ids_proxy[db_id] = {EAdapterFields.SEQ: tsvselect[0]['Coding sequence'],
                                             EAdapterFields.GENE: tsvselect[0]['Associated Gene Name'],
                                             EAdapterFields.STRAND: "-" if int(tsvselect[0]['Strand']) < 0
                                             else "+"}
        return self.ids_proxy[db_id]

    def get_variant_id_from_protein_id(self,  **kwargs):
        """
        returns all information needed to instantiate a variation

        :param trans_id: A transcript ID (either ENSAMBLE (ENS) or RefSeq (NM, XN)
        :return: list of dicts -- containing all information needed for a variant initialization
        """
        _db = kwargs.get("_db","hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")
        filter = None
        db_id = ""
        if "refseq" in kwargs:
            filter = "refseq_peptide"
            db_id = kwargs["refseq"]
        elif "ensemble" in kwargs:
            filter = "ensembl_peptide_id"
            db_id = kwargs["ensemble"]
        elif "swiss_accid" in kwargs:
            filter = "uniprot_swissprot_accession"
            db_id = kwargs["swiss_accid"]
        elif "swiss_gene" in kwargs:
            filter= "uniprot_swissprot"
            db_id = kwargs["swiss_gene"]
        else:
            warnings.war
        rq_n = self.biomart_head%(_db, _dataset) \
               + self.biomart_filter%(filter, str(db_id)) \
               + self.biomart_filter%("germ_line_variation_source", "dbSNP") \
               + self.biomart_attribute%("variation_name") \
               + self.biomart_attribute%("snp_chromosome_name") \
               + self.biomart_attribute%("chromosome_location") \
               + self.biomart_attribute%("allele") \
               + self.biomart_attribute%("snp_strand") \
               + self.biomart_attribute%("peptide_location") \
               + self.biomart_tail
        print rq_n
        tsvreader = csv.DictReader(urllib2.urlopen(self.new_biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        print tsvselect
        if not tsvselect:
            warnings.warn("No entry found for ID %s"%db_id)
            return None

        return tsvselect

    def get_protein_sequence_from_protein_id(self, **kwargs):
        """
        Returns the protein sequence for a given protein ID that can either be refeseq, uniprot or ensamble id

        :param kwargs:
        :return:
        """
        filter = None
        db_id = ""

        _db = kwargs.get("_db","hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")
        if "refseq" in kwargs:
            filter = "refseq_peptide"
            db_id = kwargs["refseq"]
        elif "ensemble" in kwargs:
            filter = "ensembl_peptide_id"
            db_id = kwargs["ensemble"]
        elif "swiss_accid" in kwargs:
            filter = "uniprot_swissprot_accession"
            db_id = kwargs["swiss_accid"]
        elif "swiss_gene" in kwargs:
            filter= "uniprot_swissprot"
            db_id = kwargs["swiss_gene"]
        else:
            warnings.warn("No correct transcript id")
            return None
        rq_n = self.biomart_head%(_db, _dataset) \
               + self.biomart_filter%(filter, str(db_id)) \
               + self.biomart_attribute%(filter) \
               + self.biomart_attribute%("peptide") \
               + self.biomart_tail
        print rq_n
        tsvreader = csv.DictReader(urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        print tsvselect
        if not tsvselect:
            warnings.warn("No entry found for ID %s"%db_id)
            return None

        return {EAdapterFields.SEQ: tsvselect[0]['Coding sequence'],
                                             EAdapterFields.GENE: tsvselect[0]['Associated Gene Name'],
                                             EAdapterFields.STRAND: "-" if int(tsvselect[0]['Strand']) < 0
                                             else "+"}

    #TODO: refactor ... function based on old code
    def get_all_variant_gene(self, locations,_db="hsapiens_gene_ensembl", _dataset='gene_ensembl_config'):
        """
        Fetches the important db ids and names for given chromosomal location
        :param chrom: integer value of the chromosome in question
        :param start: integer value of the variation start position on given chromosome
        :param stop: integer value of the variation stop position on given chromosome
        :return: The respective gene name, i.e. the first one reported
        """
        #TODO assert types
        #<!DOCTYPE Query><Query client="true" processor="TSV" limit="-1" header="1"><Dataset name="hsapiens_gene_ensembl" config="gene_ensembl_config"><Filter name="chromosomal_region" value="1:40367114:40367114,1:40702744:40702744,1:40705023:40705023,1:40771399:40771399,1:40777210:40777210,1:40881015:40881015,1:41235036:41235036,1:42048927:42048927,1:43002232:43002232,1:43308758:43308758,1:43393391:43630154,1:43772617:43772617,1:43772834:43772834" filter_list=""/><Attribute name="uniprot_genename"/></Dataset></Query>
        queries = [self.biomart_filter%("uniprot_genename", ','.join(kwargs['genes'][x:x+500])) for x in xrange(0, len(kwargs['genes']), 500)]
        # if chrom + start + stop in self.gene_proxy:
        #     return self.gene_proxy[chrom + start + stop]
        #
        # rq_n = self.biomart_head \
        #     + self.biomart_filter%("chromosome_name", str(chrom))  \
        #     + self.biomart_filter%("start", str(start))  \
        #     + self.biomart_filter%("end", str(stop))  \
        #     + self.biomart_attribute%("uniprot_genename")  \
        #     + self.biomart_tail
        #
        # tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
        # tsvselect = [x for x in tsvreader]
        # if tsvselect and tsvselect[0]:
        #     self.gene_proxy[chrom + start + stop] = tsvselect[0]['UniProt Gene Name']
        #     return tsvselect[0]['UniProt Gene Name']
        # else:
        #     logging.warn(','.join([str(chrom), str(start), str(stop)]) + ' does not denote a known gene location')
        #     return ''

    #TODO: refactor ... function based on old code
    def get_variant_ids(self, **kwargs):
        """
        Fetches the important db ids and names for given gene _or_ chromosomal location. The former is recommended.
        Result is a list of dicts with either of the tree combinations:
            - 'Ensembl Gene ID', 'Ensembl Transcript ID', 'Ensembl Protein ID'
            - 'RefSeq Protein ID [e.g. NP_001005353]', 'RefSeq mRNA [e.g. NM_001195597]', first triplet
            - 'RefSeq Predicted Protein ID [e.g. XP_001720922]', 'RefSeq mRNA predicted [e.g. XM_001125684]', first triplet
        :keyword 'chrom': integer value of the chromosome in question
        :keyword 'start': integer value of the variation start position on given chromosome
        :keyword 'stop': integer value of the variation stop position on given chromosome
        :keyword 'gene': string value of the gene of variation
        :keyword 'transcript_id': string value of the gene of variation
        :return: The list of dicts of entries with transcript and protein ids (either NM+NP or XM+XP)
        """
        # TODO type assessment
        _db = kwargs.get("_db","hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")
        ensemble_only = False
        query = None
        if len(kwargs) == 4 and 'chrom' in kwargs and 'start' in kwargs and 'stop' in kwargs and 'ensemble_only' in kwargs:
            ensemble_only = kwargs['ensemble_only']
            query = self.biomart_filter%("chromosome_name", kwargs['chrom'])  \
                + self.biomart_filter%("start", kwargs['start'])  \
                + self.biomart_filter%("end", kwargs['stop'])
        if len(kwargs) == 3 and 'chrom' in kwargs and 'start' in kwargs and 'stop' in kwargs:
            query = self.biomart_filter%("chromosome_name", kwargs['chrom'])  \
                + self.biomart_filter%("start", kwargs['start'])  \
                + self.biomart_filter%("end", kwargs['stop'])
        elif len(kwargs) == 2 and 'gene' in kwargs and 'ensemble_only' in kwargs:
            ensemble_only = kwargs['ensemble_only']
            query = self.biomart_filter%("uniprot_genename", kwargs['gene'])
        elif len(kwargs) == 2 and 'transcript_id' in kwargs and 'ensemble_only' in kwargs:
            ensemble_only = kwargs['ensemble_only']
            query = self.biomart_filter%("ensembl_transcript_id", kwargs['transcript_id'])
        elif len(kwargs) == 1 and 'gene' in kwargs:
            query = self.biomart_filter%("uniprot_genename", kwargs['gene'])
        else:
            warnings.warn("wrong arguments to get_variant_ids")

        rq_n = self.biomart_head%(_db, _dataset) \
            + query \
            + self.biomart_attribute%("ensembl_gene_id")  \
            + self.biomart_attribute%("ensembl_peptide_id")  \
            + self.biomart_attribute%("ensembl_transcript_id")  \
            + self.biomart_attribute%("strand")
        if not ensemble_only:
            rq_n += self.biomart_attribute%("refseq_mrna")  \
                + self.biomart_attribute%("refseq_peptide")
        rq_n += self.biomart_attribute%("uniprot_swissprot") + self.biomart_tail

        # logging.warning(rq_n)

        tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
        if ensemble_only:
            result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader}
        else:
            result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                  if x['RefSeq Protein ID [e.g. NP_001005353]'] and x['RefSeq mRNA [e.g. NM_001195597]']}

        if not ensemble_only:
            rq_x = self.biomart_head \
                + query  \
                + self.biomart_attribute%("ensembl_gene_id")  \
                + self.biomart_attribute%("ensembl_peptide_id")  \
                + self.biomart_attribute%("ensembl_transcript_id")  \
                + self.biomart_attribute%("refseq_peptide_predicted")  \
                + self.biomart_attribute%("refseq_mrna_predicted")  \
                + self.biomart_tail

            tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_x)).read()).splitlines(), dialect='excel-tab')

            result2 = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                       if (x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and x['RefSeq mRNA predicted [e.g. XM_001125684]'])
                       or (not x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and not x['RefSeq mRNA predicted [e.g. XM_001125684]'])}
            result.update(result2)

        g = None
        for k, v in result.iteritems():
            if 'uniprot_swissprot' in v:
                g = v['uniprot_swissprot']
        self.ids_proxy[g] = result.values()
        return result.values()

    #TODO: refactor ... function based on old code
    def get_all_variant_ids(self,  **kwargs):
        """
        Fetches the important db ids and names for given gene _or_ chromosomal location. The former is recommended.
        Result is a list of dicts with either of the tree combinations:
            - 'Ensembl Gene ID', 'Ensembl Transcript ID', 'Ensembl Protein ID'
            - 'RefSeq Protein ID [e.g. NP_001005353]', 'RefSeq mRNA [e.g. NM_001195597]', first triplet
            - 'RefSeq Predicted Protein ID [e.g. XP_001720922]', 'RefSeq mRNA predicted [e.g. XM_001125684]', first triplet
        :keyword 'locations': list of locations as triplets of integer values representing (chrom, start, stop)
        :keyword 'genes': list of genes as string value of the genes of variation
        :return: The list of dicts of entries with transcript and protein ids (either NM+NP or XM+XP)
        """
        _db = kwargs.get("_db","hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")
        end_result = dict()
        # TODO type assessment
        ensemble_only = False
        query = None
        if 'ensemble_only' in kwargs:
            ensemble_only = kwargs['ensemble_only']
        if 'locations' in kwargs:
            pass
            #TODO
            # query = self.biomart_filter%("chromosome_name", kwargs['chrom'])  \
            #         + self.biomart_filter%("start", kwargs['start'])  \
            #         + self.biomart_filter%("end", kwargs['stop'])
        elif 'genes' in kwargs:
            queries = [self.biomart_filter%("uniprot_genename", ','.join(kwargs['genes'][x:x+500])) for x in xrange(0, len(kwargs['genes']), 500)]
        else:
            logging.warning("wrong arguments to get_variant_ids")
        for query in queries:
            rq_n = self.biomart_head%(_db,_dataset) \
                + query \
                + self.biomart_attribute%("uniprot_genename")  \
                + self.biomart_attribute%("ensembl_gene_id")  \
                + self.biomart_attribute%("ensembl_peptide_id")  \
                + self.biomart_attribute%("ensembl_transcript_id")  \
                + self.biomart_attribute%("strand")
            if not ensemble_only:
                rq_n += self.biomart_attribute%("refseq_mrna")  \
                    + self.biomart_attribute%("refseq_peptide")
            rq_n += self.biomart_tail
            # rq_n += self.biomart_attribute%("uniprot_swissprot") + self.biomart_tail

            # logging.warning(rq_n)

            tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
            if ensemble_only:
                result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader}
            else:
                result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                      if x['RefSeq Protein ID [e.g. NP_001005353]'] and x['RefSeq mRNA [e.g. NM_001195597]']}
            end_result.update(result)

            if not ensemble_only:
                rq_x = self.biomart_head%(_db, _dataset) \
                    + query  \
                    + self.biomart_attribute%("uniprot_genename")  \
                    + self.biomart_attribute%("ensembl_gene_id")  \
                    + self.biomart_attribute%("ensembl_peptide_id")  \
                    + self.biomart_attribute%("ensembl_transcript_id")  \
                    + self.biomart_attribute%("refseq_peptide_predicted")  \
                    + self.biomart_attribute%("refseq_mrna_predicted")  \
                    + self.biomart_tail

                tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_x)).read()).splitlines(), dialect='excel-tab')

                for x in tsvreader:
                    if (x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and x['RefSeq mRNA predicted [e.g. XM_001125684]']):
                        end_result.setdefault(x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID'], x)
                # result2 = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                #            if (x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and x['RefSeq mRNA predicted [e.g. XM_001125684]'])
                #     this line is a troublemaker and does not help       or (not x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and not x['RefSeq mRNA predicted [e.g. XM_001125684]'])}
                # end_result.setdefault(result2)

        # g = None
        # for k, v in result.iteritems():
        #     if 'uniprot_swissprot' in v:
        #         g = v['uniprot_swissprot']
        # self.ids_proxy[g] = result.values()
        return end_result.values()
