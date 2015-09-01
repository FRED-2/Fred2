# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'mohr, walzer'

from logging import warning
from tempfile import NamedTemporaryFile

from Fred2.Core.Peptide import Peptide
from Fred2.Core.Result import Distance2SelfResult
from Fred2 import d2s


def p2s(peptides):
    if not isinstance(peptides, list):
        warning("no list of peptides given - returning None!")
        return None

    if any(not(isinstance(p, Peptide) or isinstance(p, str)) for p in peptides):
        warning("no list with peptides given - returning None!")
        return None

    revi = list()
    for x in peptides:
        if isinstance(x, Peptide):
            revi.append(str(x))
        elif isinstance(x, str):
            revi.append(x)
        else:
            warning("no list of peptides given - returning None!")
            return None
    return revi


class Distance2Self(object):
    """
        Implements calulcation routine of distance to (self) peptides
        Calculate k closest distances of peptide to peptide set represented as trie

        All our matrices have the same ordering of letters.
        If you use a new matrix, pleas make sure to use the same ordering! Otherwise the tries have to be recomputed!
    """
    def __init__(self, matrix, selfpeptides, max_res=10):
        """

        :param matrix: a DistanceMatrix
        :param selfpeptides: list of Peptides or str of Peptides all in same length
        :param max_res: max number of minimum distances returned
        """
        self.__matrix = matrix
        self.__matrixfile = NamedTemporaryFile(delete=False)
        self.__matrixfile.close()
        matrix.to_file(self.__matrixfile.name)
        self.__max_res = max_res
        #TODO @mohr: default behaviour for heterogeneous lists?
        self.__d2s = d2s.Distance2Self(self.__matrixfile.name, p2s(selfpeptides)) #TODO @mohr: what default?
        #TODO get trie from d2s and write like done with __matrixfile

    def calculate_distances(self, peptides):
        import sys
        print "@@@", sys.getrecursionlimit()
        dd = self.__d2s.getDistance2Self(p2s(peptides), self.__max_res)
        print "derp"
        result = {key: value.values() for key, value in dd.iteritems()}
        # create temporary file with peptides for distance computation

        resultDf = Distance2SelfResult.from_dict(result)
        resultDf['matrix'] = self.__matrix.name

        return resultDf