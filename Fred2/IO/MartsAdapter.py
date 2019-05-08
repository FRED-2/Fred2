# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: IO.MartsAdapter
   :synopsis: BDB-Adapter for BioMart
.. moduleauthor:: walzer, schubert
"""

import csv
import urllib.request, urllib.error, urllib.parse
import warnings
import logging
import pymysql.cursors
from operator import itemgetter

from Fred2.IO.ADBAdapter import ADBAdapter, EAdapterFields, EIdentifierTypes


class MartsAdapter(ADBAdapter):
    def __init__(self, usr=None, host=None, pwd=None, db=None, biomart=None):
        """
        Used to fetch sequences from given RefSeq id's either from BioMart if no credentials given else from a MySQLdb
co
        :param str usr: db user e.g. = 'ucsc_annot_query'
        :param str host: db host e.g. = "pride"
        :param str pwd: pw for user e.g. = 'an0q3ry'
        :param str db: db on host e.g. = "hg18_ucsc_annotation"
        """
        self.ids_proxy = dict()
        self.gene_proxy = dict()
        self.sequence_proxy = dict()

        if usr and host and pwd and db:
            self.connection = pymysql.connect(user=usr, host=host, password=pwd, db=db)
        else:
            self.connection = None

        if biomart:
            self.biomart_url = biomart
            if not self.biomart_url.endswith("/biomart/martservice?query="):
                self.biomart_url += "/biomart/martservice?query="
        else:
            self.biomart_url = "http://biomart.org/biomart/martservice?query="
            #new: http://central.biomart.org/biomart/martservice?query="
            #http://www.ensembl.org/biomart/martview/
            #http://grch37.ensembl.org/biomart/
            #http://localhost:9000/biomart/martservice?query=%3C!DOCTYPE%20Query%3E%3CQuery%20client=%22biomartclient%22%20processor=%22TSV%22%20limit=%22-1%22%20header=%221%22%3E%3CDataset%20name=%22hsapiens_gene_ensembl%22%20config=%22gene_ensembl_config%22%3E%3CFilter%20name=%22uniprot_genename%22%20value=%22TP53%22%20filter_list=%22%22/%3E%3CFilter%20name=%22germ_line_variation_source%22%20value=%22dbSNP%22%20filter_list=%22%22/%3E%3CAttribute%20name=%22snp_ensembl_gene_id%22/%3E%3CAttribute%20name=%22snp_chromosome_name%22/%3E%3CAttribute%20name=%22snp_ensembl_transcript_id%22/%3E%3CAttribute%20name=%22snp_start_position%22/%3E%3CAttribute%20name=%22snp_ensembl_peptide_id%22/%3E%3CAttribute%20name=%22snp_end_position%22/%3E%3CAttribute%20name=%22snp_external_gene_name%22/%3E%3CAttribute%20name=%22snp_strand%22/%3E%3CAttribute%20name=%22source_description%22/%3E%3CAttribute%20name=%22germ_line_variation_source%22/%3E%3CAttribute%20name=%22allele%22/%3E%3CAttribute%20name=%22variation_name%22/%3E%3C/Dataset%3E%3C/Query%3E
            #self.biomart_url = """http://134.2.9.124/biomart/martservice?query="""
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

    def get_product_sequence(self, product_id, **kwargs):
        """
        Fetches product (i.e. protein) sequence for the given id

        :param str product_id: The id to be queried
        :keyword type: Assumes given ID from type found in :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`, default is
                       ensembl_peptide_id
        :type type: :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`
        :keyword str _db: Can override MartsAdapter default db ("hsapiens_gene_ensembl")
        :keyword str _dataset: Specifies the query dbs dataset if default is not wanted ("gene_ensembl_config")

        :return: The requested sequence
        :rtype: str
        """

        _db = kwargs.get("_db", "hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")

        query_filter = "ensembl_peptide_id"
        if "type" in kwargs:
            query_id = kwargs["type"]
            if kwargs["type"] == EIdentifierTypes.REFSEQ:
                query_filter = "refseq_peptide"
            elif kwargs["type"] == EIdentifierTypes.PREDREFSEQ:
                query_filter = "refseq_peptide_predicted"
            elif kwargs["type"] == EIdentifierTypes.ENSEMBL:
                query_filter = "ensembl_peptide_id"
            else:
                logging.warn("Could not infer the origin of product id " + str(product_id))
                return None

        if product_id in self.sequence_proxy:
            return self.sequence_proxy[product_id]

        rq_n = self.biomart_head%(_db, _dataset) \
            + self.biomart_filter%(query_filter, str(product_id))  \
            + self.biomart_attribute%("peptide")  \
            + self.biomart_attribute%("external_gene_name")  \
            + self.biomart_tail

        # logging.warn(rq_n)
        tsvreader = csv.DictReader(urllib.request.urlopen(self.biomart_url
                                                   + urllib.parse.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        if not tsvselect:
            logging.warn("There seems to be no Proteinsequence for " + str(product_id))
            #print self.biomart_url+rq_n
            return None
        self.sequence_proxy[product_id] = tsvselect[0]["Peptide"][:-1] if tsvselect[0]["Peptide"].endswith('*')\
            else tsvselect[0]["Peptide"]
        return self.sequence_proxy[product_id]

    def get_transcript_sequence(self, transcript_id, **kwargs):
        """
        Fetches transcript sequence for the given id

        :param str transcript_id: The id to be queried
        :keyword type: Assumes given ID from type found in :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`, default is
                       ensembl_transcript_id
        :type type: :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`
        :keyword str _db: Can override MartsAdapter default db ("hsapiens_gene_ensembl")
        :keyword str _dataset: Specifies the query dbs dataset if default is not wanted ("gene_ensembl_config")

        :return: The requested sequence
        :rtype: str
        """

        _db = kwargs.get("_db", "hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")

        query_filter = "ensembl_transcript_id"
        if "type" in kwargs:
            query_id = kwargs["type"]
            if kwargs["type"] == EIdentifierTypes.REFSEQ:
                query_filter = "refseq_mrna"
            elif kwargs["type"] == EIdentifierTypes.PREDREFSEQ:
                query_filter = "refseq_mrna_predicted"
            elif kwargs["type"] == EIdentifierTypes.ENSEMBL:
                query_filter = "ensembl_transcript_id"
            else:
                logging.warn("Could not infer the origin of transcript id " + str(transcript_id))
                return None

        if transcript_id in self.gene_proxy:
            return self.gene_proxy[transcript_id]

        rq_n = self.biomart_head%(_db, _dataset) \
            + self.biomart_filter%(query_filter, str(transcript_id))  \
            + self.biomart_attribute%(query_filter)  \
            + self.biomart_attribute%("coding")  \
            + self.biomart_attribute%("strand")  \
            + self.biomart_tail

        tsvreader = csv.DictReader(urllib.request.urlopen(self.biomart_url +
                                                   urllib.parse.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        if not tsvselect:
            logging.warn("There seems to be no transcript sequence for " + str(transcript_id))
            #print self.biomart_url+rq_n
            return None

        self.sequence_proxy[transcript_id] = tsvselect[0]['Coding sequence']
        return self.sequence_proxy[transcript_id]

    def get_transcript_information(self, transcript_id, **kwargs):
        """
        Fetches transcript sequence, gene name and strand information for the given id

        :param str transcript_id: The id to be queried
        :keyword type: Assumes given ID from type found in :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`, default is
                       ensembl_transcript_id
        :type type: :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`
        :keyword str _db: Can override MartsAdapter default db ("hsapiens_gene_ensembl")
        :keyword str _dataset: Specifies the query dbs dataset if default is not wanted ("gene_ensembl_config")

        :return: Dictionary of the requested keys as in EAdapterFields.ENUM
        :rtype: dict
        """

        _db = kwargs.get("_db", "hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")

        query_filter = "ensembl_transcript_id"
        if "type" in kwargs:
            if kwargs["type"] == EIdentifierTypes.REFSEQ:
                query_filter = "refseq_mrna"
            elif kwargs["type"] == EIdentifierTypes.PREDREFSEQ:
                query_filter = "refseq_mrna_predicted"
            elif kwargs["type"] == EIdentifierTypes.ENSEMBL:
                query_filter = "ensembl_transcript_id"
            else:
                logging.warn("Could not infer the origin of transcript id " + str(transcript_id))
                return None

        if transcript_id in self.ids_proxy:
            return self.ids_proxy[transcript_id]

        rq_n = self.biomart_head%(_db, _dataset) \
            + self.biomart_filter%(query_filter, str(transcript_id))  \
            + self.biomart_attribute%(query_filter)  \
            + self.biomart_attribute%("coding")  \
            + self.biomart_attribute%("strand")  \
            + self.biomart_tail

        tsvreader = csv.DictReader(urllib.request.urlopen(self.biomart_url+urllib.parse.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        if not tsvselect:
            logging.warn("No Information on transcript %s"%transcript_id)
            return None

        self.ids_proxy[transcript_id] = {EAdapterFields.SEQ: tsvselect[0]['Coding sequence'],
                                                  EAdapterFields.GENE: tsvselect[0].get('Associated Gene Name', ""),
                                                  EAdapterFields.STRAND: "-" if int(tsvselect[0]['Strand']) < 0
                                                  else "+"}
        return self.ids_proxy[transcript_id]

    def get_transcript_position(self, transcript_id, start, stop, **kwargs):
        """
        If no transcript position is available for a variant, it can be retrieved if the mart has the transcripts
        connected to the CDS and the exons positions

        :param str transcript_id: The id to be queried
        :param int start: First genomic position to be mapped
        :param int stop: Last genomic position to be mapped
        :keyword type: Assumes given ID from type found in :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`, default is
                       ensembl_transcript_id
        :type type: :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`
        :keyword str _db: Can override MartsAdapter default db ("hsapiens_gene_ensembl")
        :keyword str _dataset: Specifies the query dbs dataset if default is not wanted ("gene_ensembl_config")

        :return: A tuple of the mapped positions start, stop
        :rtype: int
        """
        # ma = MartsAdapter(biomart="http://grch37.ensembl.org")
        # print ma.get_transcript_position('17953929', '17953943', 'ENST00000361221')
        # (1674, 1688)
        try:
            x = int(start)
            y = int(stop)
        except Exception as e:
            logging.warning(','.join([str(start), str(stop)]) + ' does not seem to be a genomic position.')
            return None

        _db = kwargs.get("_db", "hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")

        query_filter = "ensembl_transcript_id"
        if "type" in kwargs:
            if kwargs["type"] == EIdentifierTypes.REFSEQ:
                query_filter = "refseq_mrna"
            elif kwargs["type"] == EIdentifierTypes.PREDREFSEQ:
                query_filter = "refseq_mrna_predicted"
            elif kwargs["type"] == EIdentifierTypes.ENSEMBL:
                query_filter = "ensembl_transcript_id"
            else:
                logging.warn("Could not infer the origin of transcript id " + str(transcript_id))
                return None

        if str(start) + str(stop) + transcript_id in self.gene_proxy:
            return self.gene_proxy[str(start) + str(stop) + transcript_id]

        rq_n = self.biomart_head%(_db, _dataset) \
            + self.biomart_filter%(query_filter, str(transcript_id))  \
            + self.biomart_attribute%("exon_chrom_start")  \
            + self.biomart_attribute%("exon_chrom_end")  \
            + self.biomart_attribute%("strand")  \
            + self.biomart_attribute%("cds_start")  \
            + self.biomart_attribute%("cds_end")  \
            + self.biomart_tail

        tsvreader = csv.DictReader((urllib.request.urlopen(self.biomart_url +
                                                    urllib.parse.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
        exons = [ex for ex in tsvreader if ex["CDS start"] and ex["CDS end"]]
        cds = [dict([k, int(v)] for k, v in e.items()) for e in exons] # cast to int
        cds = sorted(cds, key=itemgetter("CDS start")) #sort by CDS Start(position in the CDS)

        cds_sum = 0
        if not cds:
            logging.warning(str(transcript_id) + ' does not seem to have exons mapped.')
            return None

        for e in cds:
            sc = e["CDS start"]
            ec = e["CDS end"]
            if not sc or not ec:
                continue
            se = e["Exon region start (bp)"]
            ee = e["Exon region end (bp)"]

            if not cds_sum < sc < ec:
                logging.warn("unable to follow the CDS, aborting genome-positional lookup in transcript!")
                print(exons)
                print(cds)
                return None
                #after sorting and filtering if this occurs points to corrupt data in mart

            if x in range(se, ee + 1):
                if not y in range(se, ee + 1):
                    logging.warning(','.join([str(start), str(stop)]) +
                                    ' spans more than one exon, aborting genome-positional lookup in transcript!')
                    return None
                else:
                    #strand dependent!!!
                    if e["Strand"] < 0:  # reverse strand!!!
                        self.gene_proxy[str(start) + str(stop) + transcript_id] =\
                            (ee - x + 1 + cds_sum, ee - y + 1 + cds_sum)
                    else:  # forward strand!!!
                        self.gene_proxy[str(start) + str(stop) + transcript_id] =\
                            (x - se + 1 + cds_sum, y - se + 1 + cds_sum)
                    return self.gene_proxy[str(start) + str(stop) + transcript_id]
            else:
                cds_sum = ec

        logging.warning(','.join([str(start), str(stop)]) + ' seems to be outside of the exons boundaries.')
        return None

    def get_gene_by_position(self, chrom, start, stop, **kwargs):
        """
        Fetches the gene name for given chromosomal location

        :param int chrom: Integer value of the chromosome in question
        :param int start: Integer value of the variation start position on given chromosome
        :param int stop: Integer value of the variation stop position on given chromosome
        :keyword str _db: Can override MartsAdapter default db ("hsapiens_gene_ensembl")
        :keyword str _dataset: Specifies the query dbs dataset if default is not wanted ("gene_ensembl_config")

        :return: The respective gene name, i.e. the first one reported
        :rtype: str
        """
        if str(chrom) + str(start) + str(stop) in self.gene_proxy:
            return self.gene_proxy[str(chrom) + str(start) + str(stop)]

        _db = kwargs.get("_db", "hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")

        rq_n = self.biomart_head%(_db, _dataset) \
            + self.biomart_filter%("chromosome_name", str(chrom))  \
            + self.biomart_filter%("start", str(start))  \
            + self.biomart_filter%("end", str(stop))  \
            + self.biomart_attribute%("external_gene_name")  \
            + self.biomart_tail

        tsvreader = csv.DictReader((urllib.request.urlopen(self.biomart_url +
                                                    urllib.parse.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        if tsvselect and tsvselect[0]:
            self.gene_proxy[str(chrom) + str(start) + str(stop)] = tsvselect[0]['Gene name']
            return self.gene_proxy[str(chrom) + str(start) + str(stop)]
        else:
            logging.warning(','.join([str(chrom), str(start), str(stop)]) + ' does not denote a known gene location')
            return ''

    def get_transcript_information_from_protein_id(self, product_id, **kwargs):
        """
        Fetches transcript sequence for the given id

        :param str product_id: The id to be queried
        :keyword type: Assumes given ID from type found in :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`, default is
                       ensembl_peptide_id
        :type type: :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`
        :keyword str _db: Can override MartsAdapter default db ("hsapiens_gene_ensembl")
        :keyword str _dataset: Specifies the query dbs dataset if default is not wanted ("gene_ensembl_config")
        :return: List of dictionary of the requested sequence, the respective strand and the associated gene name
        :rtype: list(dict)
        """

        _db = kwargs.get("_db", "hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")

        query_filter = "ensembl_peptide_id"
        if "type" in kwargs:
            if kwargs["type"] == EIdentifierTypes.REFSEQ:
                query_filter = "refseq_peptide"
            elif kwargs["type"] == EIdentifierTypes.PREDREFSEQ:
                query_filter = "refseq_peptide_predicted"
            elif kwargs["type"] == EIdentifierTypes.ENSEMBL:
                query_filter = "ensembl_peptide_id"
            else:
                logging.warn("Could not infer the origin of product id " + str(product_id))
                return None

        if product_id in self.ids_proxy:
            return self.ids_proxy[product_id]

        rq_n = self.biomart_head%(_db, _dataset) \
               + self.biomart_filter%(query_filter, str(product_id)) \
               + self.biomart_attribute%(query_filter) \
               + self.biomart_attribute%("coding") \
               + self.biomart_attribute%("external_gene_id") \
               + self.biomart_attribute%("strand") \
               + self.biomart_tail

        tsvreader = csv.DictReader(urllib.request.urlopen(self.biomart_url+urllib.parse.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        if not tsvselect:
            warnings.warn("No entry found for id %s"%product_id)
            return None

        self.ids_proxy[product_id] = {EAdapterFields.SEQ: tsvselect[0]['Coding sequence'],
                                             EAdapterFields.GENE: tsvselect[0]['Associated Gene Name'],
                                             EAdapterFields.STRAND: "-" if int(tsvselect[0]['Strand']) < 0
                                             else "+"}
        return self.ids_proxy[product_id]

    def get_variant_id_from_protein_id(self, transcript_id, **kwargs):
        """
        Returns all information needed to instantiate a variation

        :param str transcript_id: The id to be queried
        :keyword type: assumes given ID from type found in :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`, default is
                       ensembl_transcript_id
        :type type: :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`
        :keyword str _db: can override MartsAdapter default db ("hsapiens_gene_ensembl")
        :keyword str _dataset: specifies the query dbs dataset if default is not wanted ("gene_ensembl_config")

        :return: Containing all information needed for a variant initialization
        :rtype: list(dict)
        """
        _db = kwargs.get("_db", "hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")

        query_filter = "ensembl_peptide_id"
        if "type" in kwargs:
            if kwargs["type"] == EIdentifierTypes.REFSEQ:
                query_filter = "refseq_mrna"
            elif kwargs["type"] == EIdentifierTypes.PREDREFSEQ:
                query_filter = "refseq_mrna_predicted"
            elif kwargs["type"] == EIdentifierTypes.ENSEMBL:
                query_filter = "ensembl_transcript_id"
            else:
                logging.warn("Could not infer the origin of transcript id " + str(transcript_id))
                return None

        rq_n = self.biomart_head%(_db, _dataset) \
               + self.biomart_filter%(query_filter, str(transcript_id)) \
               + self.biomart_filter%("germ_line_variation_source", "dbSNP") \
               + self.biomart_attribute%("ensembl_gene_id") \
               + self.biomart_attribute%("variation_name") \
               + self.biomart_attribute%("snp_chromosome_name") \
               + self.biomart_attribute%("chromosome_location") \
               + self.biomart_attribute%("allele") \
               + self.biomart_attribute%("snp_strand") \
               + self.biomart_attribute%("peptide_location") \
               + self.biomart_tail
        #tsvreader = csv.DictReader(urllib2.urlopen(self.new_biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab') #what? brachvogel?
        tsvreader = csv.DictReader(urllib.request.urlopen(self.biomart_url+urllib.parse.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        if not tsvselect:
            warnings.warn("No entry found for id %s"%transcript_id)
            return None

        return tsvselect

    def get_ensembl_ids_from_id(self, gene_id, **kwargs):
        """
        Returns a list of gene-transcript-protein ids from some sort of id

        :param str gene_id: The id to be queried
        :keyword type: Assumes given ID from type found in list of :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes` ,
                       default is gene name
        :type type: :func:`~Fred2.IO.ADBAdapter.EIdentifierTypes`
        :keyword str _db: can override MartsAdapter default db ("hsapiens_gene_ensembl")
        :keyword str _dataset: specifies the query dbs dataset if default is not wanted ("gene_ensembl_config")

        :return: Containing information about the corresponding (linked) entries.
        :rtype: list(dict)
        """
        _db = kwargs.get("_db","hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")

        query_filter = "external_gene_name"
        if "type" in kwargs:
            if kwargs["type"] == EIdentifierTypes.HGNC:
                query_filter = "hgnc_symbol"
            elif kwargs["type"] == EIdentifierTypes.UNIPROT:
                query_filter = "uniprot_swissprot"
            elif kwargs["type"] == EIdentifierTypes.GENENAME:
                query_filter = "external_gene_name"
            elif kwargs["type"] == EIdentifierTypes.ENSEMBL:
                query_filter = "ensemble_gene_id"
            else:
                logging.warn("Could not infer the origin of gene id " + str(gene_id))
                return None

        if gene_id in self.ids_proxy:
            return self.ids_proxy[gene_id]

        rq_n = self.biomart_head%(_db, _dataset) \
               + self.biomart_filter%(query_filter, str(gene_id)) \
               + self.biomart_attribute%("ensembl_gene_id") \
               + self.biomart_attribute%("strand") \
               + self.biomart_attribute%("ensembl_transcript_id") \
               + self.biomart_attribute%("ensembl_peptide_id") \
               + self.biomart_tail
        tsvreader = csv.DictReader(urllib.request.urlopen(self.biomart_url +
                                                   urllib.parse.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        if not tsvselect:
            logging.warn("No entry found for id %s"%gene_id)
            return None

        self.ids_proxy[gene_id] = [{EAdapterFields.PROTID: gtp.get('Protein ID', ""),
                                                  EAdapterFields.GENE: gtp.get('Gene ID', ""),
                                                  EAdapterFields.TRANSID: gtp.get('Transcript ID', ""),
                                                  EAdapterFields.STRAND: "-" if int(gtp['Strand']) < 0
                                                  else "+"} for gtp in tsvselect]
        return self.ids_proxy[gene_id]

    #TODO: refactor ... function based on old code
    def get_all_variant_gene(self, locations, _db="hsapiens_gene_ensembl", _dataset='gene_ensembl_config'):
        """
        Fetches the important db ids and names for given chromosomal location

        :param int chrom: Integer value of the chromosome in question
        :param int start: Integer value of the variation start position on given chromosome
        :param int stop: Integer value of the variation stop position on given chromosome
        :return: The respective gene name, i.e. the first one reported

        """
        #TODO assert types
        #<!DOCTYPE Query><Query client="true" processor="TSV" limit="-1" header="1"><Dataset name="hsapiens_gene_ensembl" config="gene_ensembl_config"><Filter name="chromosomal_region" value="1:40367114:40367114,1:40702744:40702744,1:40705023:40705023,1:40771399:40771399,1:40777210:40777210,1:40881015:40881015,1:41235036:41235036,1:42048927:42048927,1:43002232:43002232,1:43308758:43308758,1:43393391:43630154,1:43772617:43772617,1:43772834:43772834" filter_list=""/><Attribute name="uniprot_genename"/></Dataset></Query>
        #queries = [self.biomart_filter%("uniprot_genename", ','.join(kwargs['genes'][x:x+300])) for x in xrange(0, len(kwargs['genes']), 300)]
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
        AResult is a list of dicts with either of the tree combinations:
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
            + self.biomart_attribute%("uniprot_genename")  \
            + self.biomart_attribute%("ensembl_gene_id")  \
            + self.biomart_attribute%("ensembl_peptide_id")  \
            + self.biomart_attribute%("ensembl_transcript_id")  \
            + self.biomart_attribute%("strand")
        if not ensemble_only:
            rq_n += self.biomart_attribute%("refseq_mrna")  \
                + self.biomart_attribute%("refseq_peptide")
        rq_n += self.biomart_attribute%("uniprot_swissprot") + self.biomart_tail

        # logging.warning(rq_n)

        tsvreader = csv.DictReader((urllib.request.urlopen(self.biomart_url+urllib.parse.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
        if ensemble_only:
            result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader}
        else:
            result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                  if x['RefSeq Protein ID [e.g. NP_001005353]'] and x['RefSeq mRNA [e.g. NM_001195597]']}

        if not ensemble_only:
            rq_x = self.biomart_head \
                + query  \
                + self.biomart_attribute%("uniprot_genename")  \
                + self.biomart_attribute%("ensembl_gene_id")  \
                + self.biomart_attribute%("ensembl_peptide_id")  \
                + self.biomart_attribute%("ensembl_transcript_id")  \
                + self.biomart_attribute%("refseq_peptide_predicted")  \
                + self.biomart_attribute%("refseq_mrna_predicted")  \
                + self.biomart_attribute%("strand")  \
                + self.biomart_tail

            tsvreader = csv.DictReader((urllib.request.urlopen(self.biomart_url+urllib.parse.quote(rq_x)).read()).splitlines(), dialect='excel-tab')

            result2 = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                       if (x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and x['RefSeq mRNA predicted [e.g. XM_001125684]'])
                       or (not x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and not x['RefSeq mRNA predicted [e.g. XM_001125684]'])}
            result.update(result2)

        g = None
        for k, v in result.items():
            if 'uniprot_swissprot' in v:
                g = v['uniprot_swissprot']
        self.ids_proxy[g] = list(result.values())
        return list(result.values())

    #TODO: refactor ... function based on old code
    def get_all_variant_ids(self, **kwargs):
        """
        Fetches the important db ids and names for given gene _or_ chromosomal location. The former is recommended.
        AResult is a list of dicts with either of the tree combinations:
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
            queries = [self.biomart_filter%("uniprot_genename", ','.join(kwargs['genes'][x:x+250])) for x in range(0, len(kwargs['genes']), 250)]
            logging.warning('***'+self.biomart_filter%("uniprot_genename", ','.join(kwargs['genes'][x:x+250])) for x in range(0, len(kwargs['genes']), 250))
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

            try:
                tsvreader = csv.DictReader((urllib.request.urlopen(self.biomart_url+urllib.parse.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
                if ensemble_only:
                    result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader}
                else:
                    result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                          if x['RefSeq Protein ID [e.g. NP_001005353]'] and x['RefSeq mRNA [e.g. NM_001195597]']}
                end_result.update(result)
            except:
                logging.error('Bad Mart Query: '+rq_n)

            if not ensemble_only:
                rq_x = self.biomart_head%(_db, _dataset) \
                    + query  \
                    + self.biomart_attribute%("uniprot_genename")  \
                    + self.biomart_attribute%("ensembl_gene_id")  \
                    + self.biomart_attribute%("ensembl_peptide_id")  \
                    + self.biomart_attribute%("ensembl_transcript_id")  \
                    + self.biomart_attribute%("refseq_peptide_predicted")  \
                    + self.biomart_attribute%("refseq_mrna_predicted")  \
                    + self.biomart_attribute%("strand")  \
                    + self.biomart_tail

                try:
                    tsvreader = csv.DictReader((urllib.request.urlopen(self.biomart_url+urllib.parse.quote(rq_x)).read()).splitlines(), dialect='excel-tab')

                    for x in tsvreader:
                        if (x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and x['RefSeq mRNA predicted [e.g. XM_001125684]']):
                            end_result.setdefault(x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID'], x)
                    # result2 = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                    #            if (x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and x['RefSeq mRNA predicted [e.g. XM_001125684]'])
                    #     this line is a troublemaker and does not help       or (not x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and not x['RefSeq mRNA predicted [e.g. XM_001125684]'])}
                    # end_result.setdefault(result2)
                except:
                    logging.error('Bad Mart Query: '+rq_n)

        # g = None
        # for k, v in result.iteritems():
        #     if 'uniprot_swissprot' in v:
        #         g = v['uniprot_swissprot']
        # self.ids_proxy[g] = result.values()
        return list(end_result.values())
