__author__ = 'schubert'

import unittest

from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.CleavagePrediction import CleavageSitePredictorFactory
from Fred2.EpitopeAssembly.EpitopeAssembly import EpitopeAssemblyWithSpacer


class SpacerDesignTestCase(unittest.TestCase):
    """
        Unittest for OptiTope
    """

    def setUp(self):
        epis = """GHRMAWDMM
                 VYEADDVIL""".split("\n")

        self.epis = [Peptide(x.strip()) for x in epis]
        self.alleles = [Allele("HLA-A*02:01", prob=0.5)]

    def test_standart_functions(self):
        """
        Tests default functions
        needs GLPK installed
        :return:
        """
        epi_pred = EpitopePredictorFactory("Syfpeithi")
        cl_pred = CleavageSitePredictorFactory("PCM")

        sbws = EpitopeAssemblyWithSpacer(self.epis, cl_pred, epi_pred, self.alleles, solver="cbc")
        sol = sbws.solve()
        print(sol)
        assert all(i == str(j) for i, j in zip(["GHRMAWDMM", "HH", "VYEADDVIL"], sol))

    def test_unsupported_allele_length_combination(self):
        """
        Tests default functions
        needs GLPK installed
        :return:
        """
        epi_pred = EpitopePredictorFactory("Syfpeithi")
        cl_pred = CleavageSitePredictorFactory("PCM")
        alleles = [Allele("HLA-A*02:01", prob=0.5), Allele("HLA-A*26:01", prob=0.5)]
        sbws = EpitopeAssemblyWithSpacer(self.epis, cl_pred, epi_pred, alleles, solver="cbc")
        sol = sbws.solve()
        print(sol)
        assert all(i == str(j) for i, j in zip(["GHRMAWDMM", "HH", "VYEADDVIL"], sol))

    def test_unsupported_allele_length_combination_exception(self):
        """
        Tests default functions
        needs GLPK installed
        :return:
        """
        epi_pred = EpitopePredictorFactory("Syfpeithi")
        cl_pred = CleavageSitePredictorFactory("PCM")
        alleles = [Allele("HLA-A*26:01", prob=0.5)]
        sbws = EpitopeAssemblyWithSpacer(self.epis, cl_pred, epi_pred, alleles, solver="cbc")
        self.assertRaises(ValueError, sbws.solve)


if __name__ == '__main__':
    unittest.main()
