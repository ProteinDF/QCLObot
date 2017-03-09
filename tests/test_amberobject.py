#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from qclobot.amberobject import AmberObject
from qclobot.utils import get_model
import pdfbridge

class TestAmberObject(unittest.TestCase):
    def setUp(self):
        pdb = pdfbridge.Pdb("./data/sample/4tut.addH.pdb")
        models = pdb.get_atomgroup()
        self.model = get_model(models)

    def tearDown(self):
        pass
        
    def test_opt(self):
        self.amber_obj = AmberObject("test_amber_opt")
        self.amber_obj.model = self.model
        self.amber_obj.opt()
        
    def test_md(self):
        self.amber_obj = AmberObject("test_amber_md")
        self.amber_obj.model = self.model
        self.amber_obj.md(steps=5000, dt=0.002)
        

def test_suite():
    """
    builds the test suite.
    """
    def _suite(test_class):
        return unittest.makeSuite(test_class)

    suite = unittest.TestSuite()
    suite.addTests(_suite(TestAmberObject))
    return suite

if __name__ == '__main__':
    unittest.main(defaultTest = 'test_suite')
        
    
