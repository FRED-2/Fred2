# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.Allele
   :synopsis: Allele class.
.. moduleauthor:: schubert, brachvogel, szolek, walzer

"""


from Fred2.Core.Base import MetadataLogger


class Allele(MetadataLogger):
    """
    This class represents an Allele and stores additional informations as a 
    dictionary

    :param str name: the name of the MHC allele (new nomenclature A*01:01)
    """
    def __init__(self, _name, prob=None):
        """
        :param str _name: input name in new nomenclature (A*01:01)
        """
        MetadataLogger.__init__(self)
        name = _name.split("-")[-1].replace("HLA-", "")
        self.name = name
        self.locus, rest = name.split('*')
        self.supertype, self.subtype = rest.split(':')[:2]
        self.prob = prob

        # TODO check semantics

    def __repr__(self):
        return 'HLA-%s*%s:%s' % (str(self.locus), str(self.supertype), str(self.subtype))

    def __str__(self):
        return self.name

    def __eq__(self, other):
        return self.name == other.name

    def __cmp__(self, other):
        return cmp(self.name, other.name)