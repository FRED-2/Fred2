
__author__ = 'mohr'

import time
import os
from collections import OrderedDict


class DistanceMatrix(object):

    def __init__(self, data):

        self.__matrix = data
        self.__pathToMatrixFile = None
        assert isinstance(data, OrderedDict), 'DistanceMatrix object has to be initiated with dictionary.'
        self.to_file()

    def __str__(self):

        repr = '\t'
        repr += '\t'.join(self.__matrix.keys())
        repr += '\n'
        for k in self.__matrix.keys():
            lineStr = '%s ' % k
            for kk in self.__matrix[k].keys():
                lineStr += '%s ' % self.__matrix[k][kk]

            repr += '%s\n' % lineStr
        return repr

    @property
    def path_to_matrix_file(self):

        return self.__pathToMatrixFile

    def to_file(self):

        timestr = time.strftime("%Y%m%d-%H%M%S")
        pathToMatrix = os.path.abspath("../Data/tmp/%s.mat"% timestr)
        with open(pathToMatrix,'w') as newMatrixFile:
            for k in self.__matrix.keys():
                lineStr = '%s\t' % k
                for kk in self.__matrix[k].keys():
                    lineStr += '%s\t' % self.__matrix[k][kk]
                newMatrixFile.write(lineStr + '\n')
        self.__pathToMatrixFile = pathToMatrix

    def from_file(self, filePath):

        with open(filePath, 'r') as matrixFile:
            matrixDict = OrderedDict([])
            aas = []
            distances = []
            for line in matrixFile:
                values = line.strip().split(' ')
                aas.append(values[0])
                distances.append(values[1:])

            for aa, vals in zip(aas, distances):
                matrixValues = OrderedDict([])
                for aaa, v in zip(aas,vals):
                    matrixValues[aaa] = v
                matrixDict[aa] = matrixValues
        self.__matrix = matrixDict
        self.to_file()

