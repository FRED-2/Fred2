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


class CombinedAllele(Allele):
    """
    This class represents combined HLA class II Alleles with an alpha and beta chain
    """

    def __init__(self, name, prob=None):
        """
        :param str name: input name in new nomenclature (DPA1*01:03-DPB1*01:01)
        :param float prob: optional population frequency of allele in [0,1]
        """
        MetadataLogger.__init__(self)
        alpha_chain = name.replace("HLA-", "").split("-")[0]
        beta_chain = name.replace("HLA-", "").split("-")[1]
        name = alpha_chain+"-"+beta_chain

        self.name = name
        self.alpha_chain = alpha_chain
        self.beta_chain = beta_chain
        self.alpha_locus, alpha_rest = self.alpha_chain.split('*')
        self.alpha_supertype, self.alpha_subtype = alpha_rest.split(':')[:2]
        self.beta_locus, beta_rest = self.beta_chain.split('*')
        self.beta_supertype, self.beta_subtype = beta_rest.split(':')[:2]
        self.prob = prob

    def __repr__(self):
        return 'HLA-%s*%s:%s-%s*%s:%s' % (str(self.alpha_locus), str(self.alpha_supertype), str(self.alpha_subtype),
                                          str(self.beta_locus), str(self.beta_supertype), str(self.beta_subtype))
