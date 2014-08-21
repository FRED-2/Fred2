# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'schubert', 'walzer'

import logging
#from Prediction.PSSM import Syfpeithi
from Base import MetadataLogger
from collections import OrderedDict


#class Allele(MetadataLogger):
class Allele(object):
    """
    This class represents an Allele and stores additional informations as a dictionary
    """
    def __init__(self, name):
        """
        :param name (string): the name of the MHC allele (new nomenclature A*01:01)
        """
        # MetadataLogger.__init__(self)
        self.locus, rest = name.split('*')
        self.supertype, self.subtype = rest.split(':')[:2]
        # TODO check semantics

    def __repr__(self):
        return 'HLA-%s*%s:%s' % (str(self.locus), str(self.supertype), str(self.subtype))

    def __str__(self):
        return self.__repr__()

    def to_netmhc(self, netmhc, version):
        # allele format: A0101. For netMHCpan: HLA-A01:01
        #if isinstance(netmhc, NetMHC):
            if version == 'netMHC-3.0':
                a = self.locus + self.supertype + self.subtype
                #logging.warning(a)
                if a in netmhc.mhcalleles:
                    return a
                else:
                    raise LookupError
            elif version == 'netMHCpan-2.4':
                a = 'HLA-%s%s:%s' % (self.locus, self.supertype, self.subtype)
                if a in netmhc.panalleles:
                    return a
                else:
                    raise LookupError
            else:
                raise LookupError

    def to_syfpeithi(self, syfpeithi, length):
        #if isinstance(syfpeithi, Syfpeithi):
            if '%s%s%s_%s' % (self.locus, self.supertype, self.subtype, length) not in syfpeithi.get_matrices():
                if '%s%s_%s' % (self.locus, self.supertype, length) not in syfpeithi.get_matrices():
                    raise LookupError
                else:
                    return '%s%s_%s' % (self.locus, self.supertype, length)
            else:
                return '%s%s%s_%s' % (self.locus, self.supertype, self.subtype, length)
        #else:
            #raise ValueError

