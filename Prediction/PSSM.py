# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'walzer', 'schubert'

import collections, itertools, warnings, abc
from Fred2.Core.Base import AEpitopePrediction


#DEPREcated
#
# class PSSM(dict):
#     """
#     BaseclasCATEDs for PSSM predictions.
#     """
#     def __init__(self, name, rows):
#         """
#         Initializes a PSSM prediction method object.
#         """
#         self.name = name
#         self.length = 0
#         self.max = 0
#         for row in rows:
#             self[row[0]] = map(int, row[1:])
#
#         self.length = len(row[1:])
#         #TODO: properly check lengths
#         #~ transposed = zip(*original)
#         self.max = sum(max(i) for i in map(list, zip(*self.values())))
#
#     def __str__(self):
#         return ' '.join(self.keys())  # TODO
#
#     def predict(self, peptide):
#         """
#         Make prediction for a Peptide.
#         """
#         message = "Peptide " + str(peptide.seq) + " does not match the Matrix" + self.name + " (Aminoacidcomposition or length " + str(self.length) + ")."
#         if len(peptide) != self.length or not isinstance(peptide, AASequence):
#             logging.warning(message)
#             return None
#         try:
#             score = 0
#             for i, aa in enumerate(peptide):
#                 #TODO enforce Bio.alphabet?
#                 score = score + self[aa.upper()][i]
#             return score
#         except KeyError as e:
#             logging.warn(message + ' ' + e)
#             return None
#
#     def percent_max(self, score):
#         return (100.0/float(self.max))*float(score)
#
#
# class Syfpeithi(MetadataLogger):
#
#     def __init__(self, matrix_directory='Syfpeithi'):
#         MetadataLogger.__init__(self)
#         self.matrices = dict()
#         for i in glob.glob(matrix_directory + "/*.syf"):
#             with open(i) as syf:
#                 rows = list()
#                 name = str()
#                 for line in csv.reader(syf, dialect="excel-tab"):
#                     if sum(i != '' for i in line) == 1:
#                         if len(rows) > 0 and name != '':
#                             self.matrices[name] = PSSM(name, rows)
#                         rows = list()
#                         name = str()
#                         name = line[0]
#                         #~ print '---'
#                         #~ print name
#                     elif sum(i != '' for i in line) < 1:
#                         continue
#                     else:
#                         rows.append(line)
#                         #~ print line
#
#     def predict(self, peptide, allele_name):
#         if allele_name in self.matrices:
#             return self.matrices[allele_name].predict(peptide)
#         else:
#             message = "No such allele_name: " + allele_name + "."
#             raise Exception(message)
#
#     def percent_max(self, score, allele_name):
#         if allele_name in self.matrices:
#             return self.matrices[allele_name].percent_max(score)
#         else:
#             message = "No such allele_name: " + allele_name + "."
#             raise Exception(message)
#
#     def get_matrices(self, length=None):
#         if not length:
#             return self.matrices.keys()
#         else:
#             return [i for i in self.matrices.keys() if self.matrices[i].length == length]
#
#     def make_predictions(self, peptides, alleles=None, ignore=True):
#         """
#
#         :type alleles: list of Core.Allele
#         """
#         if not alleles:
#             alleles = self.get_matrices()
#         else:
#             assert all(isinstance(a, Allele) for a in alleles), "No list of Allele"
#         assert all(isinstance(a, AASequence) for a in peptides), "No list of AASequence"
#         # TODO: refactor these 3 lines for new peptide format:
#         pepset = Core.uniquify_list(peptides, fred2_attrgetter('seq'))
#         pepset.sort(key=len)
#         pepsets = [list(g) for k, g in groupby(pepset, key=len)]
#
#         for set in pepsets:
#             for allele in alleles:
#                 try:
#                     a = allele.to_syfpeithi(self, len(set[0]))
#                 except LookupError:
#                     logging.warning("Allele not available for this method (Syfpeithi): "+str(allele))
#                     if ignore:
#                         continue
#
#                 for pepseq in set:
#                     try:
#                         syf_score = self.predict(pepseq, a)
#                         score = Score('Syfpeithi', allele, syf_score, self.percent_max(syf_score, a), None)
#                         for pep in [p for p in peptides if str(pepseq.seq) == str(p.seq)]:
#                             pep.scores.append(score)
#                     except:
#                         if not ignore:
#                             logging.warning(''.join(pepseq, " failed with ", allele))
#
#
#
#

#
#
# class BIMAS(MetadataLogger):
#     pass
#
# #~ #codesnippet for transforming old FRED hardcoded dict-list-dict matrices
#   #~ import csv
#   #~ with open('/tmp/test.tsv', 'w') as fp:
#     #~ for k in syfpeithi_matrices:
#       #~ fp.write( k+'\n')
#       #~ aas = syfpeithi_matrices[k][1].keys()
#       #~ aas.sort()
#       #~ for a in aas:
#         #~ fp.write( a +'\t')
#         #~ for i,v in enumerate(syfpeithi_matrices[k]):
#           #~ fp.write( str(syfpeithi_matrices[k][i][a]) +'\t')
#         #~ fp.write('\n')
#       #~ fp.write('\n')
#     #~ for k in syfpeithi_matrices2:
#       #~ fp.write( k+'\n')
#       #~ aas = syfpeithi_matrices2[k][1].keys()
#       #~ aas.sort()
#       #~ for a in aas:
#         #~ fp.write( a +'\t')
#         #~ for i,v in enumerate(syfpeithi_matrices2[k]):
#           #~ fp.write( str(syfpeithi_matrices2[k][i][a]) +'\t')
#         #~ fp.write('\n')
#       #~ fp.write('\n')
#

class APSSMPredictor(AEpitopePrediction):
    """
        Abstract base class for PSSM predictions.

        Implements predict functionality

    """

    def predict(self, peptides, alleles=None, **kwargs):
        """
        Returns predictions for given peptides an alleles. If no alleles are given, predictions for all available models
        are made.

        :param list(Peptide)/Peptide peptides: A single Peptide or a list of Peptides
        :param list(Alleles) alleles: a list of Alleles
        :param kwargs: optional parameter (not used yet)
        :return: Returns a Result object with the prediction results
        """
        def __load_allele_model(allele,length):
            allele_model = "%s_%s_%i"%(self.name, allele, length)
            return getattr( __import__("Fred2.Data.PSSMMatrices", fromlist=[allele_model]), allele_model)


        if isinstance(peptides, collections.Iterable):
            pep_seqs = {str(p):p for p in peptides}
        else:
            pep_seqs = {str(peptides):peptides}

        if alleles is None:
            allales_string = self.supportedAlleles
        else:
            allales_string = self.convert_alleles(alleles)

        #group peptides by length and
        result = {}
        for length, peps in itertools.groupby(pep_seqs.iterkeys(), key= lambda x: len(x)):

            #dynamicaly import prediction PSSMS for alleles and predict
            if length not in self.supportedLength:
                warnings.warn("Peptide length of %i not supported"%length, RuntimeWarning)
                continue

            for a in allales_string:
                pssm = __load_allele_model(a, length)
                ##here is the prediction and result object missing##
                for p in peps:
                    score = sum(pssm[i][p[i]] for i in xrange(length))
                    result.setdefault(str(p), []).append((p, self.name, a, score))

        return result


class Syfpeithi(APSSMPredictor):
    """
        Represents the Syfpeithi PSSM predictor
    """
    __alleles = ["A*03:01"]
    __supported_length = [9, 10, 11]
    __name = "syfpeithi"

    @property
    def name(self):
        return self.__name

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__supported_length

    def convert_alleles(self, alleles):
        pass

    def predict(self, peptides, alleles=None, **kwargs):
        #with this implementation customizations of prediction algorithm is still possible
        return super(Syfpeithi, self).predict(peptides, alleles=alleles, **kwargs)


class BIMAS(APSSMPredictor):
    """
        Represents the BIMAS PSSM predictor
    """

    __alleles = ["HLA-A*31:01"]
    __supported_length = [9, 10, 11]
    __name = "bimas"

    @property
    def name(self):
        return self.__name

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__supported_length

    def convert_alleles(self, alleles):
        return ["%s_%s%s"%(a.locus, a.supertype, a.subtype) for a in alleles]

    def predict(self, peptides, alleles=None, **kwargs):
        #with this implementation customizations of prediction algorithm is still possible
        return super(BIMAS, self).predict(peptides, alleles=alleles, **kwargs)
