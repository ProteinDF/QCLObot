#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import os
import copy
# import pprint
# import jinja2

# import proteindf_bridge as bridge

from . import __version__
from .modeler_task_edit import ModelerEdit
from .qccontrolobject import QcControlObject
from .modeler_task_protonate import QcProtonate
from .modeler_task_complement import ModelerTaskComplement
from .modeler_task_opt import ModelerTaskOpt
from .qcneutralize import QcNeutralize
from .amberobject import AmberObject
from .utils import (file2atomgroup,
                    check_format_model_list,
                    check_format_model)

import logging
logger = logging.getLogger(__name__)


class QcModeler(QcControlObject):
    def __init__(self):
        super(QcModeler, self).__init__()
        self.global_tasks = {}

    def show_version(self):
        logger.info("=" * 80)
        logger.info("QcModeler version: {version}".format(
            version=str(__version__)))
        logger.info("=" * 80)

    def _run_task_cmd(self, task):
        task = copy.deepcopy(task)
        # pprint.pprint(task)
        # name = task.pop("name")
        task_name = task.get("name")
        # print("name:", name)

        if "edit" in task.keys():
            # pprint.pprint(task['edit'])
            self.global_tasks[task_name] = self._run_edit(task)
        elif "protonate" in task.keys():
            self.global_tasks[task_name] = self._run_protonate(task)
        elif "neutralize" in task.keys():
            self.global_tasks[task_name] = self._run_neutralize(task)
        elif "complement" in task.keys():
            self.global_tasks[task_name] = self._run_complement(task)
        elif "opt" in task.keys():
            self.global_tasks[task_name] = self._run_opt(task)
        elif "md" in task.keys():
            self.global_tasks[task_name] = self._run_md(task_name, task['md'])
        else:
            logger.warning("not found task command: {cmds}".format(
                cmds=str([x for x in task.keys()])))

    def _run_edit(self, task):
        logger.info("run edit")
        modeler_edit = ModelerEdit(self, task)
        modeler_edit.run()
        modeler_edit.finalize()

        return modeler_edit

    def _run_protonate(self, task):
        logger.info("run protonate")

        modeler_protonate = QcProtonate(self, task)
        modeler_protonate.run()
        modeler_protonate.finalize()

        return modeler_protonate

    def _run_neutralize(self, task):
        logger.info("run neutralize")

        modeler_neutralize = QcNeutralize(self, task)
        modeler_neutralize.run()
        modeler_neutralize.finalize()

        return modeler_neutralize

    def _run_complement(self, task):
        logger.info("run complement")

        modeler_task_complement = ModelerTaskComplement(self, task)
        modeler_task_complement.run()
        modeler_task_complement.finalize()

        return modeler_task_complement

    def _run_opt(self, task):
        logger.info("run opt")

        modeler_task_opt = ModelerTaskOpt(self, task)
        modeler_task_opt.run()
        modeler_task_opt.finalize()

        return modeler_task_opt

    def _run_md(self, name, args):
        logger.info("run md")
        assert(isinstance(name, str))
        assert(isinstance(args, dict))

        steps = 1
        if "steps" in args:
            steps = args["steps"]
        dt = 0.002
        if "dt" in args:
            dt = args["dt"]

        amber = self._run_amber_setup(name, args)
        amber.md(steps=steps,
                 dt=dt)

        dest = args.get("dest", None)
        if dest:
            amber.write_output_model(dest)

        return amber

    def _run_amber_setup(self, name, args):
        amber = AmberObject(name=name)
        model = self._get_input_model(args)
        amber.model = model

        if "solvation" in args:
            solvation_args = args["solvation"]

            amber.solvation_method = "cap"
            if "method" in solvation_args:
                amber.solvation_method = solvation_args["method"]

            if "model" in solvation_args:
                amber.solvation_model = solvation_args["model"]

        if "belly_mask" in args:
            amber.use_belly = True
            for bellymask_target in args["belly_mask"]:
                bellymask_target = bellymask_target.lower()
                if bellymask_target == "water":
                    amber.bellymask_WAT = True
                if bellymask_target == "ions":
                    amber.bellymask_ions = True

        return amber

    def _get_input_model(self, args):
        """入力model(AtomGroup)を返す

        'reference'が指定されている場合は、そのtaskオブジェクトの出力(AtomGroup)を、
        'src'が指定されている場合は、そのファイルから作成されたAtomGroupを返す
        """
        answer = None
        if "reference" in args:
            ref_task = self.get_reference_task(args["reference"])
            answer = ref_task.output_model
        elif "src" in args:
            src = args.get("src", None)
            if src is None:
                logger.critical('not found "src"')
                raise
            atomgroup = file2atomgroup(src)
            if check_format_model_list(atomgroup):
                logger.info("this model is MODELS. pickup first model.")
                models = atomgroup
                if models.get_number_of_groups() > 0:
                    model_list = [k for k, v in models.groups()]
                    assert(len(model_list) > 0)
                    answer = models.get_group(model_list[0])
            else:
                answer = atomgroup
        else:
            logger.critical(
                "NOT found input model. use \"reference\" or \"src\" command.")
        assert(check_format_model(answer))

        return answer

    def get_reference_task(self, reference_task_name):
        assert(isinstance(reference_task_name, str))
        if reference_task_name in self.global_tasks.keys():
            return self.global_tasks[reference_task_name]
        else:
            logger.critical("unknown task name: {name}".format(
                name=reference_task_name))
            raise
