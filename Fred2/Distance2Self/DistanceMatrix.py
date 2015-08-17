
__author__ = 'mohr'

import datetime
import os
from collections import OrderedDict


class DistanceMatrix(object):

    def __init__(self, data, saveMatrixFile=False):

        self.__saveMatrixFile = saveMatrixFile
        self.__matrix = data
        self.__pathToMatrixFile = None
        assert isinstance(data, OrderedDict), 'DistanceMatrix object has to be initiated with dictionary.'
        self.to_file()

    def __del__(self):

        if not self.__saveMatrixFile:
            os.remove(self.__pathToMatrixFile)

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

    @property
    def path_to_matrix_file(self):
        return self.__pathToMatrixFile

    def to_file(self):
        timestr = datetime.datetime.now().strftime("%y%m%d_%H%M%S.%f")
        current = os.path.join(os.path.dirname(__file__))
        pathToMatrix = os.path.join(current,'..',"Data/tmp/%s.mat" % timestr)

        with open(pathToMatrix,'w') as newMatrixFile:
            for k in self.__matrix.keys():
                lineStr = '%s ' % k
                for kk in self.__matrix[k].keys():
                    lineStr += '%s ' % self.__matrix[k][kk]
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

