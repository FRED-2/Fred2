# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.Allele
   :synopsis: HLA Allele class.
.. moduleauthor:: schubert, brachvogel, szolek, walzer

"""


from Fred2.Core.Base import MetadataLogger


class Allele(MetadataLogger):
    """
    This class represents an HLA Allele and stores additional information
    """

    def __init__(self, name, prob=None):
        """
        :param str name: input name in new nomenclature (A*01:01)
        :param float prob: optional population frequency of allele in [0,1]
        """
        MetadataLogger.__init__(self)
        name = name.split("-")[-1].replace("HLA-", "")
        self.name = name
        self.locus, rest = name.split('*')
        self.supertype, self.subtype = rest.split(':')[:2]
        self.prob = prob

    def __repr__(self):
        return 'HLA-%s*%s:%s' % (str(self.locus), str(self.supertype), str(self.subtype))

    def __str__(self):
        return self.name

    def __eq__(self, other):
        return str(self.name) == str(other)

    def __cmp__(self, other):
        return cmp(self.name, str(other))