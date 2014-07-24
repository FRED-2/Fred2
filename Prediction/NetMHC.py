# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'szolek', 'walzer'

import time
import os
import warnings

from operator import attrgetter

from Core.Base import MetadataLogger, AASequence, Score
from Core.Allele import Allele
from Core import Utils


class NetMHC(MetadataLogger):
    def __init__(self, matrix_directory='Syfpeithi'):
        MetadataLogger.__init__(self)
        self.netmhc_path = None
        self.netpan_path = None

    def make_predictions(self, peptides, alleles=None, method='netMHC-3.0', ignore=True):
        if not alleles:
            alleles = self.get_matrices()
        else:
            assert all(isinstance(a, Allele) for a in alleles), "No list of Allele"
            alleles = [a.to_netmhc() for a in alleles]
        assert all(isinstance(a, AASequence) for a in peptides), "No list of AASequence"
        pepset = Utils.uniquify_list(peptides, attrgetter('seq'))

        runid = time.time()
        peptide_infile = os.path.join(tempdir, 'peptides_%s.txt' % runid)
        prediction_outfile = os.path.join(tempdir, 'predictions_%s.txt' % runid)

        write_peptide_file(pepset, peptide_infile)

        #MW: I would suggest a netmhc predictor object having a default path, alleles and doing the errorhandling?
        os.system(netmhc_path + ' -a %s -p %s > %s' % (','.join(allele), peptide_infile, prediction_outfile))

        with open(prediction_outfile, 'r') as g:
            i_separator = 1
            for line in g:
                if line.startswith('-------------------'):
                    i_separator += 1
                    continue

                if i_separator%3 == 0:  # rows between the second and 3rd ----- separator lines
                    sep = line.split()
                    pepseq = None
                    if ('WB' in sep or 'SB' in sep) and len(sep)==7 :
                        pepseq=sep[1]
                        score_triplet = ('NetMHC', sep[6], sep[2])
                        affinity_triplet = ('NetMHC', sep[6], sep[3])
                        rank_triplet = ('NetMHC', sep[6], sep[4])
                    elif len(sep)==6:
                        pepseq=sep[1]
                        score_triplet = ('NetMHC', sep[5], sep[2])
                        affinity_triplet = ('NetMHC', sep[5], sep[3])
                        rank_triplet = ('NetMHC', sep[5], 'NB')
                    else:
                        if not ignore:
                            print "netMHC interpretation failed with ", line
                        continue
                    #~ _, pepseq, score, affinity = line.split()[:4]
                    #~ score_triplet = ('netMHC', allele, float(score))
                    #~ affinity_triplet = ('netMHC', allele, float(affinity))

                    for peptide in pepset[pepseq]:
                        peptide.log_metadata('score', score_triplet)
                        peptide.log_metadata('affinity', affinity_triplet)


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
