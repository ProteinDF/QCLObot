#!/usr/bin/env python
# -*- coding: utf-8 -*-

import yaml
import logging
logger = logging.getLogger(__name__)
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import proteindf_bridge as bridge


class QcControlObject(object):
    def __init__(self):
        self._senarios = []
        self._vars = {}

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
        for d in yaml.load_all(contents, Loader=yaml.SafeLoader):
            self._senarios.append(d)

    def _save_yaml(self, data, path):
        assert(isinstance(data, dict))
        with open(path, 'w') as f:
            yaml.dump(data, f, encoding='utf8', allow_unicode=True)


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
        logger.info("QcControlObject::show_version()")
