# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'schubert', 'walzer'

from Core.Base import MetadataLogger


class Allele(MetadataLogger):
    """
    This class represents an Allele and stores additional informations as a dictionary
    """
    def __init__(self, name):
        """
        :param name: (String) the name of the MHC allele (new nomenclature A*01:01)
        """
        MetadataLogger.__init__(self)
        # MetadataLogger.__init__(self)
        name = name.split("-")[-1]
        self.name = name
        self.locus, rest = name.split('*')
        self.supertype, self.subtype = rest.split(':')[:2]

        # TODO check semantics

    def __repr__(self):
        return 'HLA-%s*%s:%s' % (str(self.locus), str(self.supertype), str(self.subtype))