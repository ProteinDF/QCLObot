#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pprint

import proteindf_bridge as bridge

from .taskobject import TaskObject
from .utils import check_format_model_list, check_format_model, file2atomgroup

import logging
logger = logging.getLogger(__name__)


class ModelerTaskObject(TaskObject):
    """Task Object for Modeler
    """

    def __init__(self, parent, task):
        from .qcmodeler import QcModeler
        assert(isinstance(parent, QcModeler))
        assert(isinstance(task, dict))

        self._parent = parent

        super(ModelerTaskObject, self).__init__(**task)

        # print(">>>> ModelerTaskObject::__init__()")
        # print("self._data:")
        # pprint.pprint(self._data)
        # print("parent:")
        # pprint.pprint(parent)

        self._data.update(task)
        # print("self._data(update):")
        # pprint.pprint(self._data)
        # print("<<<<")
        self.model = self._get_input_model(self._data)

    def __del__(self):
        self._write_output_model()

    # -------------------------------------------------------------------------
    # model (input atomgroup)
    #   "model" data is an AtomGroup object formatted by 'MODEL', not 'protein'.

    def _get_model(self):
        answer = None
        model_raw_data = self._data.get("model", None)
        if model_raw_data is not None:
            answer = bridge.AtomGroup(model_raw_data)
        return answer

    def _set_model(self, model):
        if check_format_model(model):
            self._data['model'] = model.get_raw_data()
        else:
            logger.critical("not support the format; use model format")
            raise
    model = property(_get_model, _set_model)

    # --------------------

    def _get_input_model(self, args):
        """入力model(AtomGroup)を返す

        'reference'が指定されている場合は、そのtaskオブジェクトの出力(AtomGroup)を、
        'src'が指定されている場合は、そのファイルから作成されたAtomGroupを返す
        """
        answer = None
        # print(">>>> _get_input_model()")
        # pprint.pprint(args)
        if "reference" in args:
            ref_task = self._parent.get_reference_task(args["reference"])
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
        # pprint.pprint(answer)
        assert(check_format_model(answer))

        # print("<<<< _get_input_model()")
        return answer

    def _get_reference_task(self, reference_task_name):
        assert(isinstance(reference_task_name, str))
        if reference_task_name in self.global_tasks.keys():
            return self.global_tasks[reference_task_name]
        else:
            logger.critical("unknown task name: {name}".format(
                name=reference_task_name))
            raise

    def _write_output_model(self):
        output_path = self._data.get("dest", None)

        if output_path is not None:
            self._atomgroup2file(self.output_model, output_path)

    def _atomgroup2file(self, atomgroup, output_path):
        '''output_pathの拡張子に応じて、bridge形式またはpdb形式でatomgroupをファイルに書き出す
        '''
        assert(isinstance(atomgroup, bridge.AtomGroup))

        abspath = os.path.abspath(output_path)
        (basename, ext) = os.path.splitext(abspath)
        ext = ext.lower()
        if ext in (".pdb", ".ent"):
            self._atomgroup2pdb(atomgroup, abspath)
        else:
            logger.info("save {path} as bridge file.".format(path=abspath))
            bridge.save_msgpack(atomgroup.get_raw_data(), abspath)

    def _atomgroup2pdb(self, atomgroup, pdbfile,
                       model_name="model_1"):
        '''atomgroupをpdb形式で出力する

        atomgroupがmodels(複数のmodelで構成されている)の場合はそのまま出力する。
        atomgroupがmodelの場合は、model_nameを付加して出力する。
        '''
        assert(isinstance(pdbfile, str))

        pdb = bridge.Pdb(mode='amber')

        protein = atomgroup
        if check_format_model(atomgroup):
            # transform MODEL object to the protein(models)
            # which has only one model.
            protein = bridge.AtomGroup()
            protein.set_group(model_name, atomgroup)

        pdb.set_by_atomgroup(protein)
        with open(pdbfile, 'w') as f:
            f.write(str(pdb))

    # output_model -----------------------------------------------------
    def _get_output_model(self):
        answer = None
        model_raw_data = self._data.get("output_model", None)
        if model_raw_data is not None:
            answer = bridge.AtomGroup(model_raw_data)
        return answer

    def _set_output_model(self, model):
        assert(check_format_model(model))
        self._data['output_model'] = model.get_raw_data()
    output_model = property(_get_output_model, _set_output_model)
