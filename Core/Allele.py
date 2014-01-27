'''
Created on Apr 24, 2013

These classes represent allele and allele sets, and are equally design as the peptide and peptideSet class

@author: Benjamin Schubert
'''

from trunk.Fred2.model.Base import MetadataLogger
from collections import OrderedDict


class Allele(MetadataLogger):
    '''
    This class represents an Allele and stores additional informations as a dictionary
    '''


    def __init__(self, name):
        '''
        Constructor
        
        @param name (string): the name of the MHC allele (new nomenclature A*01:01)
        '''
        MetadataLogger.__init__(self)
        self.name = name
        
        
    def __repr__(self):
        return self.name
    
    
class AlleleSet(MetadataLogger, OrderedDict):
    '''
        This class stores several alleles as a dictionary key by allele name!
    '''

    def __init__(self, alleles=None):
        MetadataLogger.__init__(self)
        OrderedDict.__init__(self)
        if alleles is not None:
            for allele in alleles:
                self[allele.name]=allele


    def filter(self, filter_criterion):
        '''
            function for filtering the alleleSet
            @param filter_function (Func): A function accepting two values. First is the name of the allele, second is the allele object
                @type filter_function: Is a anonymous function returning true or false
            @return: a list of allele objects full filling the filter_criterion
                @rtype: list
        '''
        return filter(filter_criterion, self.items())
        

                
                
            