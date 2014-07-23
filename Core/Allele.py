# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'schubert', 'walzer'

import warnings
from Prediction import PSSM
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
        :param name (string): the name of the MHC allele (new nomenclature A*01:01)
        """
        MetadataLogger.__init__(self)
        self.locus, rest = name.split('*')
        self.supertype, self.subtype = rest.split(':')[:2]
        # TODO check semantics

    def __repr__(self):
        return 'HLA-%s%s:%s' % (self.locus, self.supertype, self.subtype)

    def __str__(self):
        self.__repr__()

    def to_netmhc(self, version):
        if version == 'netmhc':
            return self.locus + self.supertype + self.subtype
        elif version == 'netmhcpan':
            return 'HLA-%s%s:%s' % (self.locus, self.supertype, self.subtype)
        else:
            raise ValueError
        # TODO
        #warnings.warn("HLA-ID not known to method " + hla_id)
        #    net = ['A0101', 'A0201', 'A0202', 'A0203', 'A0204', 'A0206', 'A0211', 'A0212', 'A0216', 'A0219', 'A0301', 'A1101',
        #   'A2301', 'A2402', 'A2403', 'A2601', 'A2602', 'A2902', 'A3001', 'A3002', 'A3101', 'A3301', 'A6801', 'A6802',
        #   'A6901', 'B0702', 'B0801', 'B0802', 'B1501', 'B1801', 'B2705', 'B3501', 'B3901', 'B4001', 'B4002', 'B4402',
        #   'B4403', 'B4501', 'B5101', 'B5301', 'B5401', 'B5701', 'B5801']

    def to_syfpeithi(self, syfpeithi, length):
        if isinstance(syfpeithi, PSSM.Syfpeithi):
            if '%s%s%s_%s' % (self.locus, self.supertype, self.subtype, length) not in syfpeithi.get_matrices():
                if '%s%s_%s' % (self.locus, self.supertype, length) not in syfpeithi.get_matrices():
                    raise LookupError
                else:
                    return '%s%s_%s' % (self.locus, self.supertype, length)
            else:
                return '%s%s%s_%s' % (self.locus, self.supertype, self.subtype, length)
        else:
            raise ValueError


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
