# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
from tempfile import NamedTemporaryFile

__author__ = 'mohr'

import os
import subprocess
import time
import warnings


from Fred2.Core.Result import Distance2SelfResult
from Fred2.Data import DistanceMatrices
from Fred2.Distance2Self.DistanceMatrix import DistanceMatrix


class Distance2Self(object):
    """
        Implements calulcation routine of distance to (self) peptides
        Calculate k closest distances of peptide to peptide set represented as trie

        All our matrices have the same ordering of letters.
        If you use a new matrix, pleas make sure to use the same ordering! Otherwise the tries have to be recomputed!
    """


    def __init__(self, _matrix, trie=None):

        self.__matrix = _matrix
        self.__trie = trie

        # Where should files reside ?
        self.__externalPathDistanceCalculator = '/path/to/compute_distances_ivac'
        self.__externalPathTrieGenerator = '/path/to/get_TrieArray'


    def generate_trie(self, peptidesFile, outfile='peptideTrie', peptideLength=9):

        cmd = self.__externalPathTrieGenerator + " %s %s %s %s"

        timestr = time.strftime("%Y%m%d-%H%M%S")
        pathToTrie = os.path.abspath("../Data/tries/%s_%s.trie"% (outfile, timestr))

        self.__trie = pathToTrie
        subprocess.check_output(cmd%(peptidesFile, self.__matrix.path_to_matrix_file, peptideLength, pathToTrie), shell=True)

    def calculate_distances(self, peptides, pathToTrie=None, n=10):

        if pathToTrie is None:
            trie = self.__trie
        else:
            trie = pathToTrie
            warnings.warn(
                "Order of amino acids in distance matrix which has been used to construct trie has to be the same in "
                "distance computation!",UserWarning)

        # create temporary file with peptides for distance computation
        tmpFile = NamedTemporaryFile(delete=False)

        with open(tmpFile.name, "w") as peptidesFile:
            for pep in peptides:
                peptidesFile.write('%s\n' % pep)

        cmd = self.__externalPathDistanceCalculator + " %s %s %s %s"
        result = self.parse_external_result(
            subprocess.check_output(cmd % (self.__matrix.path_to_matrix_file, trie, tmpFile.name, n), shell=True))
        os.remove(tmpFile.name)

        return result

    def parse_external_result(self, result):

        """

        :rtype : DataFrame
        """
        parsedResult = {}

        for line in result.strip().split('\n'):
            splitted = line.strip().split(" ")[-1].split(";")
            distanceValues = []
            peptide = splitted[0].split(":")[0]

            for s in splitted[:-1]:
                distanceValues.append(float(s.split(",")[-1]))

            parsedResult[peptide] = distanceValues

        resultDf = Distance2SelfResult.from_dict(parsedResult)
        resultDf['trie'] = self.__trie.split('/')[-1]
        resultDf['matrix'] = self.__matrix.path_to_matrix_file.split('/')[-1]

        return resultDf
