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
        self.header = ["pos", "peptide", "logscore", "affinity(nM)", "Bind Level", "Protein Name", "Allele"]

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
                a = allele.to_netmhc()
            except ValueError:
                logging.warn("Allele not available for netMHC")
                continue
            cmd = self.netmhc_path + ' -a %s -p %s > %s' % (a, tmp_file)
            result = subprocess.check_output(cmd, shell=True)

            netsplit = [x.lstrip().split() for x in result.split('\n')[11:-3]]
            result = dict()
            for i in netsplit:
                if len(i) == len(self.header):
                    result[i[1]] = dict(zip(self.header, i))
                else:
                    result[i[1]] = dict(zip(self.header[:4]+self.header[5:], i))

            for p in peptides:
                if p.seq in result:
                    p.scores.append(Score('netMHC-3.0', allele, result[p.seq]['logscore'], result[p.seq]['affinity(nM)'], None))




    def netmhcpan_predictions(pepset, allele, netmhc_path, tempdir=Configuration.tempdir):
        # allele format: A0101. For netMHCpan: HLA-A01:01

        allele = convert_hla_id(allele, 'netmhcpan')

        runid = time.time()
        peptide_infile = os.path.join(tempdir, 'peptides_%s.txt' % runid)
        prediction_outfile = os.path.join(tempdir, 'predictions_%s.txt' % runid)

        write_peptide_file(pepset, peptide_infile)

        os.system(netmhc_path + ' -a %s -p %s > %s' % (allele, peptide_infile, prediction_outfile))

        with open(prediction_outfile, 'r') as g:
            i_separator = 0
            for line in g:
                if line.startswith('-------------------'):
                    i_separator += 1
                    continue

                if i_separator == 2:  # rows between the second and 3rd ----- separator lines
                    _, _, pepseq, _, score, affinity, rank = line.split()[:7]
                    score_triplet = ('netMHCpan', allele, float(score))
                    affinity_triplet = ('netMHCpan', allele, float(affinity))
                    rank_triplet = ('netMHCpan', allele, float(rank))

                    for peptide in pepset[pepseq]:
                        peptide.log_metadata('score', score_triplet)
                        peptide.log_metadata('affinity', affinity_triplet)
                        peptide.log_metadata('rank', rank_triplet)
