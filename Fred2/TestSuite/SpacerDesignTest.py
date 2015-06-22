__author__ = 'schubert'

import unittest


from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.CleavagePrediction import CleavageSitePredictorFactory
from Fred2.EpitopeAssembly.EpitopeAssembly import EpitopeAssemblyWithSpacer

class OptiTopeTestCase(unittest.TestCase):
    """
        Unittest for OptiTope
    """

    def setUp(self):
        epis ="""GHRMAWDMM
                 VYEADDVIL""".split("\n")

        self.epis = map(lambda x: Peptide(x.strip()),epis)
        self.alleles =[Allele("HLA-A*02:01",prob=0.5)]

    def testStandartFunctions(self):
        """
        Tests default functions
        needs GLPK installed
        :return:
        """
        epi_pred =  EpitopePredictorFactory("Syfpeithi")
        cl_pred = CleavageSitePredictorFactory("PCM")

        sbws = EpitopeAssemblyWithSpacer(self.epis,cl_pred,epi_pred,self.alleles)
        sol = sbws.approximate()
        print sol
        assert all(i == str(j) for i,j in zip(["GHRMAWDMM","HH","VYEADDVIL"],sol))

if __name__ == '__main__':
    unittest.main()