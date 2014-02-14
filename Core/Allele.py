"""
Created on Apr 24, 2013

These classes represent allele and allele sets, and are equally design as the peptide and peptideSet class

@author: Benjamin Schubert
"""
from Fred2.Prediction import PSSM
from Base import MetadataLogger
from collections import OrderedDict


def convert_hla_id(hla_id, format, length=9):
    try:
        locus, rest = hla_id.split('*')
        supertype, subtype = rest.split(':')[:2]

        if isinstance(format, PSSM.Syfpeithi):
            if '%s%s%s_%s' % (locus, supertype, subtype, length) not in format.get_matrices():
                if '%s%s_%s' % (locus, supertype, length) not in format.get_matrices():
                    raise LookupError
                else:
                    return '%s%s_%s' % (locus, supertype, length)
            else:
                return '%s%s%s_%s' % (locus, supertype, subtype, length)
        if format == 'netmhc':
            return locus + supertype + subtype
        elif format == 'netmhcpan':
            return 'HLA-%s%s:%s' % (locus, supertype, subtype)
        else:
            raise ValueError

    except ValueError as e:  # if format is unrecognized or hla-ID can't be parsed, just return what we've got
        warnings.warn("Could not recognize HLA-ID " + hla_id)
        return hla_id  #why return just what we got?! it won't work anyway. Hm, let's think about it
    except LookupError as e:
        warnings.warn("HLA-ID not known to method " + hla_id)
        return '%s_%s%s_%s' % (locus, supertype, subtype, length)


class Allele(MetadataLogger):
    """
    This class represents an Allele and stores additional informations as a dictionary
    """

    def __init__(self, name):
        """
        Constructor
        
        @param name (string): the name of the MHC allele (new nomenclature A*01:01)
        """
        MetadataLogger.__init__(self)
        self.name = name

    def __repr__(self):
        return self.name


class AlleleSet(MetadataLogger, OrderedDict):
    """
        This class stores several alleles as a dictionary key by allele name!
    """

    def __init__(self, alleles=None):
        MetadataLogger.__init__(self)
        OrderedDict.__init__(self)
        if alleles is not None:
            for allele in alleles:
                self[allele.name] = allele

    def filter(self, filter_criterion):
        """
            function for filtering the alleleSet
            @param filter_function (Func): A function accepting two values. First is the name of the allele, second is the allele object
                @type filter_function: Is a anonymous function returning true or false
            @return: a list of allele objects full filling the filter_criterion
                @rtype: list
        """
        return filter(filter_criterion, self.items())
