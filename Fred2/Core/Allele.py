# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Core.Allele
   :synopsis: HLA Allele class.
.. moduleauthor:: schubert, brachvogel, szolek, walzer

"""


from Fred2.Core.Base import MetadataLogger


class AlleleFactory(type):
    def __call__(cls, name, prob=None):
        if cls is Allele:
            if name.startswith("H2-"):
                return MouseAllele(name, prob=prob)
            if name.count("*") > 1:
                return CombinedAllele(name, prob=prob)
        return type.__call__(cls, name, prob=prob)


class Allele(MetadataLogger):
    """
    This class represents an HLA Allele and stores additional information
    """
    __metaclass__ = AlleleFactory

    def __init__(self, name, prob=None):
        """
        :param str name: input name in new nomenclature (A*01:01)
        :param float prob: optional population frequency of allele in [0,1]
        """
        MetadataLogger.__init__(self)
        name = name.split("-")[-1].replace("HLA-", "")
        self.organism = "HLA"
        self.name = name
        self.locus, rest = name.split('*')
        self.supertype, self.subtype = rest.split(':')[:2]
        self.prob = prob

    def __repr__(self):
        return '%s-%s*%s:%s' % (str(self.organism), str(self.locus), str(self.supertype), str(self.subtype))

    def __str__(self):
        return repr(self)

    def __eq__(self, other):
        return str(self) == str(other)

    def __cmp__(self, other):
        return cmp(str(self), str(other))

    def __hash__(self):
        return hash(repr(self))


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
        self.organism = "HLA"
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

    @property
    def locus(self):
        return self.alpha_locus, self.beta_locus

    @property
    def supertype(self):
        return self.alpha_supertype, self.beta_supertype

    @property
    def subtype(self):
        return self.alpha_subtype, self.beta_subtype


class MouseAllele(Allele):
    """
    This class represents a mouse MHC allele with the following nomenclature:
    H2-Xxx

    http://www.imgt.org/IMGTrepertoireMHC/LocusGenes/index.php?repertoire=listIG_TR&species=mouse&group=MHC
    """
    def __init__(self, name, prob=None):
        MetadataLogger.__init__(self)
        allele = name.split("-")[-1].replace("H-2-", "")

        self.organism = "H-2"
        self.name = allele
        self.prob = prob

    def __repr__(self):
        return '%s-%s%s%s' % (str(self.organism), str(self.locus), str(self.supertype), str(self.subtype))

    @property
    def locus(self):
        return self.name[0]

    @property
    def supertype(self):
        return self.name[1]

    @property
    def subtype(self):
        return self.name[2:]
