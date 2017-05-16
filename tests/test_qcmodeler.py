#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import cProfile
from pstats import Stats

from qclobot.qcmodeler import QcModeler
import pdfbridge

use_profiler = False

class TestQcModeler(unittest.TestCase):
    def setUp(self):
        if use_profiler:
            self.pr = cProfile.Profile()
            self.pr.enable()

    def tearDown(self):
        if use_profiler:
            p = Stats (self.pr)
            p.strip_dirs()
            p.sort_stats ('cumtime')
            p.print_stats()
        
    def test_run(self):
        modeler = QcModeler()
        modeler.run("./data/sample/modeler.yml")

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
    print("check: RUN")
    unittest.main(defaultTest = 'test_suite')
        
    
