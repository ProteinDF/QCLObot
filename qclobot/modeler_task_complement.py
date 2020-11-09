#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .modeler_task_object import ModelerTaskObject
from .amberobject import AmberObject
# from .qcerror import QcTaskError

import logging
logger = logging.getLogger(__name__)


class ModelerTaskComplement(ModelerTaskObject):
    _task_name = "complement"

    def __init__(self, parent, task):
        super().__init__(parent, task)

    def run(self):
        self.cd_workdir()

        amber = AmberObject(model=self.model, work_dir=self.work_dir)
        amber.complement()
        self.output_model = amber.output_model

        self.restore_cwd()
        return 0
