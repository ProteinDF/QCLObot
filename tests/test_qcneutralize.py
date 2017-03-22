#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from qclobot.qcneutralize import QcNeutralize
import pdfbridge

class TestNeutralize(unittest.TestCase):
    def setUp(self):
        pdb = pdfbridge.Pdb("./data/sample/GDG_H.pdb")
        models = pdb.get_atomgroup()
        self.model = models["model_1"]
        
    def tearDown(self):
        pass

    def test_neutralize(self):
        neutralize = QcNeutralize(name="test_neutralize")
        neutralize.model = self.model
        
        retcode = neutralize.run()
        self.assertEqual(retcode, 0)

        neutralize.atomgroup2pdb(neutralize.output_model,
                                 "GDG_x.pdb")

def test_suite():
    """
    builds the test suite.
    """
    def _suite(test_class):
        return unittest.makeSuite(test_class)

    suite = unittest.TestSuite()
    suite.addTests(_suite(TestNeutralize))
    return suite


if __name__ == '__main__':
    unittest.main(defaultTest = 'test_suite')
        
