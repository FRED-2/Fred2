
__author__ = 'mohr'

import datetime
import os
from collections import OrderedDict


class DistanceMatrix(object):

    def __init__(self, data, name="DistanceMatrix"):
        self.__matrix = data
        self.name = name
        assert isinstance(data, OrderedDict), 'DistanceMatrix object has to be initiated with dictionary.'

    def __str__(self):
        repr = '\t'
        repr += '\t'.join(self.__matrix.keys())
        repr += '\n'
        for k in self.__matrix.keys():
            lineStr = '%s\t' % k
            for kk in self.__matrix[k].keys():
                lineStr += '%s\t' % self.__matrix[k][kk]

            repr += '%s\n' % lineStr
        return repr

    def to_file(self, file_path):
        with open(file_path, 'w') as newMatrixFile:
            for k in self.__matrix.keys():
                lineStr = '%s ' % k
                for kk in self.__matrix[k].keys():
                    lineStr += '%s ' % self.__matrix[k][kk]
                newMatrixFile.write(lineStr + '\n')

    def from_file(self, file_path):
        with open(file_path, 'r') as matrixFile:
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

