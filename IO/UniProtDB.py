# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'walzer'

import csv
import re
import urllib2
import warnings
import logging

from Bio import SwissProt


class UniProtDB():
    def load(self, file, sprot="sprot"):
        if sprot != "trembl":
            sprot = "sprot"

        try:
            sprot_handle = open(file, 'r')
        except:
            nuf = self.file_pre + sprot + "_" + file + self.file_suf
            try:
                sprot_handle = open(nuf, 'r')
            except:
                es = "file or species " + file + " not found! Pelase check that the uniprot file exists!"
                raise IOError(es)

        self.sprot_map = [[record.gene_name, record.sequence] for record in SwissProt.parse(sprot_handle)]

    def __init__(self, dbdir=None):
        self.file_suf = '.dat'
        if dbdir:
            self.file_pre = str(dbdir) + '/uniprot_'
        else:
            self.file_pre = '/share/databases/uniprot/knowledgebase/taxonomic_divisions/uniprot_'

        self.sprot_map = []

    def getGeneSprot(self, gene):
        result = list()
        if not self.sprot_map:
            warnings.warn("db empty! Pelase check that the uniprot file gets loaded!")
        for rec in self.sprot_map:
            if gene in rec[0]:
                result.append(rec[1])
        return result

