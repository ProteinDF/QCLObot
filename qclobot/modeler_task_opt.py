#!/usr/bin/env python
# -*- coding: utf-8 -*-


from .modeler_task_object import ModelerTaskObject
from .amberobject import AmberObject

import logging
logger = logging.getLogger(__name__)


class ModelerTaskOpt(ModelerTaskObject):
    ''' execute opt
    '''

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
        self.cd_workdir("Opt")

        self.engine.model = self.model

        opt = self._data.get('opt')
        if "solvation" in opt:
            solvation_args = opt["solvation"]

            self.engine.solvation_method = "cap"
            if "method" in solvation_args:
                self.engine.solvation_method = solvation_args["method"]

            if "model" in solvation_args:
                self.engine.solvation_model = solvation_args["model"]

        if "restraint" in opt:
            self.engine.use_restraint = True
            if "weight" in opt["restraint"]:
                self.engine.restraint_weight = opt["restraint"]["weight"]
            if "mask" in opt["restraint"]:
                self.engine.restraint_mask = opt["restraint"]["mask"]

        if "belly_mask" in opt:
            self.engine.use_belly = True
            for bellymask_target in opt["belly_mask"]:
                bellymask_target = bellymask_target.lower()
                if bellymask_target == "water":
                    self.engine.bellymask_WAT = True
                if bellymask_target == "ions":
                    self.engine.bellymask_ions = True
                if bellymask_target == "H":
                    self.engine.bellymask_H = True

        self.engine.opt()

        self.output_model = self.engine.output_model
        self.restore_cwd()

        return self
