#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import copy
import pprint
import yaml

try:
    import msgpack
except:
    import msgpack_pure as msgpack
import jinja2

import proteindf_bridge as bridge

from ._version import __version__
from .qcprotonate import QcProtonate
from .qcneutralize import QcNeutralize
from .amberobject import AmberObject
from .utils import (file2atomgroup,
                    check_format_models,
                    check_format_model)

import logging
logger = logging.getLogger(__name__)


class QcControl_Base(object):
    def __init__(self):
        self._senarios = []
        self._vars = []
        self.global_tasks = {}

    def run(self, input_path):
        self.show_version()
        self._load_yaml(input_path)

        # exec senarios
        for senario in self._senarios:
            if 'vars' in senario:
                self._run_vars(senario['vars'])
            if 'tasks' in senario:
                tasks = senario['tasks']
                for task in tasks:
                    self._run_task(task)
            else:
                logger.warn('NOT found "tasks" section')


    def _load_yaml(self, input_path):
        with open(input_path) as f:
            contents = f.read()
        contents = bridge.Utils.to_unicode(contents)

        self._senarios = []
        for d in yaml.load_all(contents):
            self._senarios.append(d)


    # ------------------------------------------------------------------
    # vars
    # ------------------------------------------------------------------
    def _run_vars(self, in_vars_data):
        assert(isinstance(in_vars_data, dict))

        self._vars = dict(in_vars_data)

    # ------------------------------------------------------------------
    # task
    # ------------------------------------------------------------------
    def _run_task(self, task):
        assert(isinstance(task, dict))

        if "name" not in task:
            logger.critical("name is not defined.")
        task_name = task.get("name")
        logger.info("--- TASK: {} ---".format(task_name))

        # setup default value
        # self._set_default_to_task(self._tasks['default'], task)

        # exec command
        self._run_task_cmd(task)

    def _set_default_to_task(self, ref_dic, out_dic):
        pass

    def _run_task_cmd(self, task):
        """execute task command

        override this function.
        """
        return None

    # ------------------------------------------------------------------
    # others
    # ------------------------------------------------------------------
    def show_version(self):
        """show version
        """
        pass
        logger.info("QcControl_Base::show_version()")


# ----------------------------------------------------------------------
class QcModeler(QcControl_Base):
    def __init__(self):
        super(QcModeler, self).__init__()

    def show_version(self):
        logger.info("=" * 80)
        logger.info("QcModeler version: {version}".format(version=str(__version__)))
        logger.info("=" * 80)

    def _run_task_cmd(self, task):
        task = copy.deepcopy(task)
        name = task.pop("name")

        if "protonate" in task.keys():
            self.global_tasks[name] = self._run_protonate(name, task['protonate'])
        elif "ff" in task.keys():
            self.global_tasks[name] = self._run_ff(name, task['ff'])
        elif "opt" in task.keys():
            self.global_tasks[name] = self._run_opt(name, task['opt'])
        elif "md" in task.keys():
            self.global_tasks[name] = self._run_md(name, task['md'])
        elif "neutralize" in task.keys():
            self.global_tasks[name] = self._run_neutralize(name, task['neutralize'])
        else:
            logger.warning("not found task command: {cmds}".format(cmds=str([ x for x in task.keys()])))


    def _run_protonate(self, name, args):
        logger.info("run prptonate")
        assert(isinstance(name, str))
        assert(isinstance(args, dict))

        backend = "reduce"
        prot_obj = QcProtonate(name=name,
                               backend=backend)

        input_model = self._get_input_model(args)
        prot_obj.model = input_model

        prot_obj.run()
        assert(check_format_model(prot_obj.output_model))

        dest = args.get("dest", None)
        if dest:
            prot_obj.write_output_model(dest)

        return prot_obj


    def _run_neutralize(self, name, args):
        logger.info("run neutralize")
        assert(isinstance(name, str))
        assert(isinstance(args, dict))

        neutralize_obj = QcNeutralize(name=name)

        input_model = self._get_input_model(args)
        neutralize_obj.model = input_model

        neutralize_obj.run()
        assert(check_format_model(neutralize_obj.output_model))

        dest = args.get("dest", None)
        if dest:
            neutralize_obj.write_output_model(dest)

        return neutralize_obj


    def _run_ff(self, name, args):
        logger.info("check ff")
        assert(isinstance(name, str))
        assert(isinstance(args, dict))

        amber = self._run_amber_setup(name, args)
        amber.check_force_field()

        return amber

    def _run_opt(self, name, args):
        logger.info("run opt")
        assert(isinstance(name, str))
        assert(isinstance(args, dict))

        amber = self._run_amber_setup(name, args)
        amber.opt()

        dest = args.get("dest", None)
        if dest:
            amber.write_output_model(dest)

        return amber


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
        amber.md(steps = steps,
                 dt = dt)

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
            ref_task = self._get_reference_task(args["reference"])
            answer = ref_task.output_model
        elif "src" in args:
            src = args.get("src", None)
            if src == None:
                logger.critical('not found "src"')
                raise
            atomgroup = file2atomgroup(src)
            if check_format_models(atomgroup):
                logger.info("this model is MODELS. pickup first model.")
                models = atomgroup
                if models.get_number_of_groups() > 0:
                    model_list = [ k for k, v in models.groups() ]
                    assert(len(model_list) > 0)
                    answer = models.get_group(model_list[0])
            else:
                answer = atomgroup
        else:
            logger.critical("NOT found input model. use \"reference\" or \"src\" command.")
        assert(check_format_model(answer))

        return answer


    def _get_reference_task(self, reference_task_name):
        assert(isinstance(reference_task_name, str))
        if reference_task_name in self.global_tasks.keys():
            return self.global_tasks[reference_task_name]
        else:
            logger.critical("unknown task name: {name}".format(name = reference_task_name))
            raise
