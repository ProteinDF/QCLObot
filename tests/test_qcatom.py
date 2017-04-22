#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from qclobot.qcatom import QcAtom

class TestQcAtom(unittest.TestCase):
    def setUp(self):
        pass
        
    def tearDown(self):
        pass
                                                

    def test_DZVP(self):
        H = QcAtom("H")
        H.basisset = "O-DZVP.H"
        C = QcAtom("C")
        C.basisset = "O-DZVP.C"

        self.assertEqual(H.basisset, "O-DZVP.H")
        self.assertEqual(H.get_number_of_AOs(), 2)
        self.assertEqual(C.basisset, "O-DZVP.C")
        self.assertEqual(C.get_number_of_AOs(), 14)

    def test_DZVP2(self):
        H = QcAtom("H")
        H.basisset = "O-DZVP2.H"
        C = QcAtom("C")
        C.basisset = "O-DZVP2.C"

        self.assertEqual(H.basisset, "O-DZVP2.H")
        self.assertEqual(H.get_number_of_AOs(), 5)
        self.assertEqual(C.basisset, "O-DZVP2.C")
        self.assertEqual(C.get_number_of_AOs(), 14)
        
def test_suite():
    """
    builds the test suite.
    """
    def _suite(test_class):
        return unittest.makeSuite(test_class)
    
    suite = unittest.TestSuite()
    suite.addTests(_suite(TestQcAtom))
    return suite


if __name__ == '__main__':
    unittest.main(defaultTest = 'test_suite')
    
