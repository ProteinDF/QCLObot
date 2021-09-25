#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import proteindf_bridge as bridge

from .taskobject import TaskObject

import logging
logger = logging.getLogger(__name__)


class QcRemoveWAT(TaskObject):
    ''' remove WAT(HOH) groups
    '''

    def __init__(self, name):
        """ initialize QcRemoveWAT object
        """
        # initialize base object
        super().__init__(name=name)

    ####################################################################
    # method
    ####################################################################
    def run(self, output_path=""):
        self.cd_workdir()

        model = bridge.Utils.remove_WAT(self.model)
        self.output_model = model

        self.restore_cwd()
        return True
