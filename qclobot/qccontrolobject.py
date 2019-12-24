#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import yaml
import pprint
import logging
logger = logging.getLogger(__name__)
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import proteindf_bridge as bridge

from .qcerror import QcError, QcControlError

class QcControlObject(object):
    def __init__(self):
        self._senarios = []
        self._senarios_index = 0

        self._vars = {}
        self._tasks = []
        self._tasks_index = 0

    def run(self, input_path):
        self.show_version()
        self._load_yaml(input_path)

        # exec senarios
        self._senarios_index = 0
        while self._senarios_index < len(self._senarios):
            senario = self._senarios[self._senarios_index]
            if 'vars' in senario:
                self._run_vars(senario['vars'])
            if 'tasks' in senario:
                self._tasks.extend(senario['tasks'])
                self._run_tasks()
            else:
                logger.warn('NOT found "tasks" section')

            self._senarios_index += 1

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
        #logger.info("VARS []".format())
        self._vars.update(in_vars_data)

    # ------------------------------------------------------------------
    # task
    # ------------------------------------------------------------------
    def _run_tasks(self):
        self._tasks_index = 0
        while self._tasks_index < len(self._tasks):
            logger.debug("[run_tasks] begin")
            logger.debug("[run_tasks] {}/{}".format(self._tasks_index, len(self._tasks)))
            logger.debug(pprint.pformat(self._tasks, width=2))
            logger.debug("[run_tasks] end")

            task = self._tasks[self._tasks_index]
            #logger.info("[RUN TASK] {}/{}".format(self._task_index, len(self._tasks)))
            #logger.info(pprint.pformat(task))

            self._run_task(task)

            self._tasks_index += 1
            # if retval == None:
            # else:
            #     if isinstance(retval, dict):
            #         if "include_tasks" in retval:
            #             tasks = retval.get("include_tasks")
            #             self._tasks[task_index:task_index] = copy.deepcopy(tasks)
            #     else:
            #         raise(QcError("program error"))

    def _run_task(self, task):
        """execute task

        """
        if not isinstance(task, dict):
            raise QcControlError("task type mismatch", task)

        if "name" not in task:
            logger.critical("name is not defined.")
        task_name = task.get("name")
        logger.info("TASK [{}]".format(task_name))

        # setup default value
        # self._set_default_to_task(self._tasks['default'], task)

        # exec command
        answer = self._run_task_cmd(task)
        return answer

    def _set_default_to_task(self, ref_dic, out_dic):
        pass

    def _run_task_cmd(self, task):
        """execute task command

        override this function.
        """
        return None

    def _exec_task_debug(self, task):
        assert(isinstance(task, dict))

        is_break = False
        if "debug" in task:
            logger.info("==== debug ====")
            contents = task.get("debug", {})
            msg = contents.get("msg", None)
            if msg:
                msg = pprint.pformat(msg)
                logger.info(msg)

            logger.info("===============")

        return is_break

    def _exec_include_tasks(self, in_task):
        assert(isinstance(in_task, dict))

        re_enter = False
        include_tasks = in_task.get("include_tasks", None)
        if include_tasks != None:
            #in_task.pop("include_tasks")
            if isinstance(include_tasks, str):
                include_path = str(include_tasks)
                logger.info("#include tasks from: {}".format(include_path))

                # read YAML file
                insert_tasks = []
                with open(include_path) as f:
                    contents = f.read()
                    contents = bridge.Utils.to_unicode(contents)

                    for d in yaml.load_all(contents, Loader=yaml.SafeLoader):
                        insert_tasks.extend(d)

                # insert tasks
                self._tasks.pop(self._tasks_index)
                #self._tasks.insert(self._tasks_index, insert_tasks)
                self._tasks[self._tasks_index : self._tasks_index] = insert_tasks
                
                # re-enter task
                self._tasks_index -= 1

                logger.debug("[include_tasks] begin insert of the following contents")
                logger.debug("self._tasks_index: {}".format(self._tasks_index))
                logger.debug(pprint.pformat(self._tasks))
                logger.debug("[include_tasks] end of insert")

                re_enter = True

            else:
                raise QcControlError("the value of include_tasks was not string.", include_tasks)

        return re_enter

    # ------------------------------------------------------------------
    # others
    # ------------------------------------------------------------------
    def show_version(self):
        """show version
        """
        pass
        logger.info("QcControlObject::show_version()")
