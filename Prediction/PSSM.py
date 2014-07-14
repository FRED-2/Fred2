# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
import re
import warnings
import glob
import csv

from Core.Peptide import Peptide, PeptideSet
from Core.Base import MetadataLogger

from Bio.Seq import Seq
from Bio.Alphabet import generic_protein



class PSSM(dict):
  """
  Baseclass for PSSM predictions.
  """

  def __init__(self,name,rows):
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
    self.max = sum(max(i) for i in map(list,zip(*self.values())))

  def predict(self,Peptide):
    """
    Make prediction for a Peptide.
    """
    message = "Peptide does not match the Matrix" + self.name + " (Aminoacidcomposition or length " + str(self.length) + ")."
    if len(Peptide) != self.length:
      warnings.warn(message)
      return None
    try:
      score = 0
      for i,aa in enumerate(str(Peptide)):
        #TODO enforce Bio.alphabet
        score = score + self[aa.upper()][i]
      return score
    except KeyError:
      warnings.warn(message)
      return None
    return score
    
  def percent_max(self, score):
    return (100.0/float(self.max))*float(score)
    

class Syfpeithi(MetadataLogger):

  def __init__(self, matrix_directory = 'Syfpeithi'):
    MetadataLogger.__init__(self)
    self.matrices = dict()
    for i in glob.glob(matrix_directory + "/*.syf"):
      with open(i) as syf:
        rows = list()
        name = str()
        for line in csv.reader(syf, dialect="excel-tab"):
          if sum(i != '' for i in line) == 1:
            if len(rows) > 0 and name != '':
              self.matrices[name] = PSSM(name,rows)
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
      message =  "No such allele_name: " + allele_name + "."
      raise Exception(message)
      return None

  def percent_max(self, score, allele_name):
    if allele_name in self.matrices:
      return self.matrices[allele_name].percent_max(score)
    else:
      message =  "No such allele_name: " + allele_name + "."
      raise Exception(message)
      return None

  def get_matrices(self, length = None):
    if not length:
      return self.matrices.keys()
    else:
      return [i for i in self.matrices.keys() if self.matrices[i].length==length]


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
  