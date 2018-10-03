#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import cProfile
from pstats import Stats

from qclobot.amberobject import AmberObject
from qclobot.utils import get_model
import proteindf_bridge as bridge

class TestAmberObject(unittest.TestCase):
    def setUp(self):
        pass
        #self.pr = cProfile.Profile()
        #self.pr.enable()

        pdb = bridge.Pdb("./data/sample/4tut.addH.pdb")
        models = pdb.get_atomgroup()
        self.model = get_model(models)

    def tearDown(self):
        #p = Stats (self.pr)
        #p.strip_dirs()
        #p.sort_stats ('cumtime')
        #p.print_stats()
        pass

    def test_opt(self):
        self.amber_obj = AmberObject("test_amber_opt")
        self.amber_obj.model = self.model
        self.amber_obj.opt()

    def test_md(self):
        self.amber_obj = AmberObject("test_amber_md")
        self.amber_obj.model = self.model
        self.amber_obj.md(steps=100, dt=0.002)


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
