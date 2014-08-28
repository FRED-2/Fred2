# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'walzer'

import logging
import glob
import csv
from operator import attrgetter
from itertools import groupby

from Core.Base import MetadataLogger, AASequence, Score, fred2_attrgetter
from Core.Allele import Allele
import Core


class PSSM(dict):
    """
    Baseclass for PSSM predictions.
    """
    def __init__(self, name, rows):
        """
        Initializes a PSSM prediction method object.
        """
        self.name = name
        self.length = 0
        self.max = 0
        for row in rows:
            self[row[0]] = map(int, row[1:])

        self.length = len(row[1:])
        #TODO: properly check lengths
        #~ transposed = zip(*original)
        self.max = sum(max(i) for i in map(list, zip(*self.values())))

    def __str__(self):
        return ' '.join(self.keys())  # TODO

    def predict(self, peptide):
        """
        Make prediction for a Peptide.
        """
        message = "Peptide " + str(peptide.seq) + " does not match the Matrix" + self.name + " (Aminoacidcomposition or length " + str(self.length) + ")."
        if len(peptide) != self.length or not isinstance(peptide, AASequence):
            logging.warning(message)
            return None
        try:
            score = 0
            for i, aa in enumerate(peptide):
                #TODO enforce Bio.alphabet?
                score = score + self[aa.upper()][i]
            return score
        except KeyError as e:
            logging.warn(message + ' ' + e)
            return None

    def percent_max(self, score):
        return (100.0/float(self.max))*float(score)


class Syfpeithi(MetadataLogger):

    def __init__(self, matrix_directory='Syfpeithi'):
        MetadataLogger.__init__(self)
        self.matrices = dict()
        for i in glob.glob(matrix_directory + "/*.syf"):
            with open(i) as syf:
                rows = list()
                name = str()
                for line in csv.reader(syf, dialect="excel-tab"):
                    if sum(i != '' for i in line) == 1:
                        if len(rows) > 0 and name != '':
                            self.matrices[name] = PSSM(name, rows)
                        rows = list()
                        name = str()
                        name = line[0]
                        #~ print '---'
                        #~ print name
                    elif sum(i != '' for i in line) < 1:
                        continue
                    else:
                        rows.append(line)
                        #~ print line

    def predict(self, peptide, allele_name):
        if allele_name in self.matrices:
            return self.matrices[allele_name].predict(peptide)
        else:
            message = "No such allele_name: " + allele_name + "."
            raise Exception(message)

    def percent_max(self, score, allele_name):
        if allele_name in self.matrices:
            return self.matrices[allele_name].percent_max(score)
        else:
            message = "No such allele_name: " + allele_name + "."
            raise Exception(message)

    def get_matrices(self, length=None):
        if not length:
            return self.matrices.keys()
        else:
            return [i for i in self.matrices.keys() if self.matrices[i].length == length]

    def make_predictions(self, peptides, alleles=None, ignore=True):
        """

        :type alleles: list of Core.Allele
        """
        if not alleles:
            alleles = self.get_matrices()
        else:
            assert all(isinstance(a, Allele) for a in alleles), "No list of Allele"
        assert all(isinstance(a, AASequence) for a in peptides), "No list of AASequence"
        # TODO: refactor these 3 lines for new peptide format:
        pepset = Core.uniquify_list(peptides, fred2_attrgetter('seq'))
        pepset.sort(key=len)
        pepsets = [list(g) for k, g in groupby(pepset, key=len)]

        for set in pepsets:
            for allele in alleles:
                try:
                    a = allele.to_syfpeithi(self, len(set[0]))
                except LookupError:
                    logging.warning("Allele not available for this method (Syfpeithi): "+str(allele))
                    if ignore:
                        continue

                for pepseq in set:
                    try:
                        syf_score = self.predict(pepseq, a)
                        score = Score('Syfpeithi', allele, syf_score, self.percent_max(syf_score, a), None)
                        for pep in [p for p in peptides if str(pepseq.seq) == str(p.seq)]:
                            pep.scores.append(score)
                    except:
                        if not ignore:
                            logging.warning(''.join(pepseq, " failed with ", allele))


class BIMAS(MetadataLogger):
    pass

#~ #codesnippet for transforming old FRED hardcoded dict-list-dict matrices
  #~ import csv
  #~ with open('/tmp/test.tsv', 'w') as fp:
    #~ for k in syfpeithi_matrices:
      #~ fp.write( k+'\n')
      #~ aas = syfpeithi_matrices[k][1].keys()
      #~ aas.sort()
      #~ for a in aas:
        #~ fp.write( a +'\t')
        #~ for i,v in enumerate(syfpeithi_matrices[k]):
          #~ fp.write( str(syfpeithi_matrices[k][i][a]) +'\t')
        #~ fp.write('\n')
      #~ fp.write('\n')
    #~ for k in syfpeithi_matrices2:
      #~ fp.write( k+'\n')
      #~ aas = syfpeithi_matrices2[k][1].keys()
      #~ aas.sort()
      #~ for a in aas:
        #~ fp.write( a +'\t')
        #~ for i,v in enumerate(syfpeithi_matrices2[k]):
          #~ fp.write( str(syfpeithi_matrices2[k][i][a]) +'\t')
        #~ fp.write('\n')
      #~ fp.write('\n')
  