# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'szolek', 'walzer'

import time
import os
import warnings

import Prediction

from IO.IO import convert_hla_id
from Core.Peptide import Peptide, PeptideSet
#import Fred2.toolbox.Configuration as Configuration  -----  Discuss what to do with configuration files! 


def netmhc_predictions(pepset, alleles = None, netmhc_path = 'netMHC/netMHC-3.0', tempdir=Configuration.tempdir, ignore = True):
    # allele format: A0101. For netMHCpan: HLA-A01:01

    allele = list()
    #~ allele = convert_hla_id(allele, 'netmhc')
    #~ allele in netmhc3
    net = ['A0101','A0201','A0202','A0203','A0204','A0206','A0211','A0212','A0216','A0219','A0301','A1101','A2301','A2402','A2403','A2601','A2602','A2902','A3001','A3002','A3101','A3301','A6801','A6802','A6901','B0702','B0801','B0802','B1501','B1801','B2705','B3501','B3901','B4001','B4002','B4402','B4403','B4501','B5101','B5301','B5401','B5701','B5801']

    if not alleles:
        allele = net
    else:
        allele = [convert_hla_id(i, 'netmhc') for i in alleles]
    #~ for a in alleles:
        #~ if not ignore and str(a) not in alla: ignore switch see syfpeithi function for example
        #~ print "netMHC failed with ", str(a) #todo should be warning
        #~ allele.append(str(a))

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


def syfpeithi_predictions(pepset, alleles = None, ignore = True):
    syf = Prediction.Syfpeithi()
    if not alleles:
        alleles = syf.get_matrices()
    else:
        alleles = [convert_hla_id(i, syf) for i in alleles]
    for pepseq in pepset:
        for a in alleles:
            try:
                affinity = syf.predict(pepseq,a)
                affinity_triplet = ('Syfpeithi', a, affinity)
                score = syf.predict(pepseq,a)
                score_triplet = ('Syfpeithi', a, syf.percent_max(score,a))
                for pep in pepset[pepseq]:
                    pep.log_metadata('affinity', affinity_triplet)
                    pep.log_metadata('score', score_triplet)
            except:
                if not ignore:
                    warnings.warn(''.join(pepseq, " failed with ", a))
