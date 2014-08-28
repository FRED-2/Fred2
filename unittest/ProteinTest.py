"""Unit test for protein class
"""

import unittest
from Fred2.Core import Protein

class TestProteinClass(unittest.TestCase):
    def setUp(self):
        # generate a Protein to test it
        self.protein = Protein.Protein("ASDERWQTGHKILPMNVFCY", 'someID')

    def test_get_gene_id(self):
        self.assertEqual(self.protein.gene_id, "someID",
                         'incorrect Id')

    def test_get_string(self):
        self.assertEqual(str(self.protein), "ASDERWQTGHKILPMNVFCY",
                         'incorrect sequence')

if __name__ == '__main__':
    unittest.main()
