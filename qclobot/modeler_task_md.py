#!/usr/bin/env python
# -*- coding: utf-8 -*-


from .modeler_task_object import ModelerTaskObject
from .amberobject import AmberObject

import logging
logger = logging.getLogger(__name__)


class ModelerTaskMd(ModelerTaskObject):
    ''' execute MD
    '''
    _task_name = "md"

    def __init__(self, parent, task):
        super().__init__(parent, task)
        self._engine = 'amber'
        self._engine_obj = None

    def _get_engine(self):
        if self._engine_obj is None:
            if self._engine == 'amber':
                self._engine_obj = AmberObject(model=self.model, work_dir=self.work_dir)

        return self._engine_obj
    engine = property(_get_engine)

    def run(self):
        self.cd_workdir("MD")

        self.engine.model = self.model

        md = self._data.get('md')

        steps = 1
        if "steps" in md:
            steps = md["steps"]
        dt = 0.002
        if "dt" in md:
            dt = md["dt"]

        self.engine.md(steps=steps, dt=dt)

        self.output_model = self.engine.output_model
        self.restore_cwd()

        return self
