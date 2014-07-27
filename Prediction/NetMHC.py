# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'szolek', 'walzer'

import re
import logging
import subprocess

from tempfile import NamedTemporaryFile

from Core.Base import MetadataLogger, AASequence, Score
from Core.Allele import Allele
import Core
import IO


class NetMHC(MetadataLogger):
    def __init__(self, netmhc=None, netpan=None):
        MetadataLogger.__init__(self)
        self.netmhc_path = netmhc
        self.netpan_path = netpan
        self.mhcheader = ["pos", "peptide", "logscore", "affinity(nM)", "Bind Level", "Protein Name", "Allele"]
        self.panheader = ["pos", "Allele", "peptide", "Protein Name", "logscore", "affinity(nM)", "%Rank", "BindLevel"]

    def make_predictions(self, peptides, alleles=None, method='netMHC-3.0', ignore=True):
        if not alleles:
            return
        else:
            assert all(isinstance(a, Allele) for a in alleles), "No list of Allele"
        assert all(isinstance(a, AASequence) for a in peptides), "No list of AASequence"
        pepset = Core.uniquify_list(peptides, Core.fred2_attrgetter('seq'))
        tmp_file = NamedTemporaryFile(delete=True)
        IO.write_peptide_file(pepset, tmp_file)

        for allele in alleles:
            try:
                a = allele.to_netmhc(method)
            except ValueError:
                logging.warn("Allele not available for netMHC")
                continue
            if method == 'netMHC-3.0':
                cmd = self.netmhc_path + ' -a %s %s' % (a, tmp_file)
            else:
                cmd = self.netpan_path + ' -a %s -p %s > %s' % (a, tmp_file)
            result = subprocess.check_output(cmd, shell=True)

            netsplit = [x.lstrip().split() for x in result.split('\n')[11:-3]] if method == 'netMHC-3.0' \
                else [x.lstrip().split() for x in result.split('\n')[55:-6]]
            result = dict()
            if method == 'netMHC-3.0':
                for i in netsplit:
                    if len(i) == len(self.mhcheader):
                        result[i[1]] = dict(zip(self.mhcheader, i))
                    else:
                        result[i[1]] = dict(zip(self.mhcheader[:4]+self.mhcheader[5:], i))
            else:
                for i in netsplit:
                        result[i[1]] = dict(zip(self.panheader, i))  # here no matter if i is shorter than header

            for p in peptides:
                if p.seq in result:
                    p.scores.append(Score(method, allele, result[p.seq]['logscore'], result[p.seq]['affinity(nM)'], None))
