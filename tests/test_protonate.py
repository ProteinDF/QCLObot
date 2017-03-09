#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from qclobot.qcprotonate import QcProtonate
from qclobot.utils import get_model
import pdfbridge

class TestProtonate(unittest.TestCase):
    def setUp(self):
        pdb = pdfbridge.Pdb("./data/sample/4tut.pdb")
        self.models = pdb.get_atomgroup()
        self.protonate = QcProtonate(name="test_protonate",
                                     atomgroup=self.models)
        
    def tearDown(self):
        pass

    def test_protonate(self):
        retcode = self.protonate.protonate()
        self.assertEqual(retcode, 0)

        self.protonate.protonate_group()
        

def test_suite():
    """
    builds the test suite.
    """
    def _suite(test_class):
        return unittest.makeSuite(test_class)

    suite = unittest.TestSuite()
    suite.addTests(_suite(TestProtonate))
    return suite


if __name__ == '__main__':
    unittest.main(defaultTest = 'test_suite')
        
