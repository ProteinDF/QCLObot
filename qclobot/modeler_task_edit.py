#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import copy
# import pprint

import proteindf_bridge as bridge

from .modeler_task_object import ModelerTaskObject
# from .qcerror import QcTaskError

import logging
logger = logging.getLogger(__name__)


class ModelerEdit(ModelerTaskObject):
    def __init__(self, parent, task):
        super(ModelerEdit, self).__init__(parent, task)

        assert('edit' in task.keys())
        # self._data['edit'] = copy.deepcopy(task['edit'])

    # ==================================================================
    # method
    # ==================================================================
    def run(self):
        # logger.info(">>>> {0[0]}:{0[1]}:{0[2]}".format(bridge.locate()))
        new_model = bridge.AtomGroup(self.model)
        # pprint.pprint(self._data)
        # print(">>>> model")
        # print(self.model)

        edit = self._data.get('edit')
        keys = edit.keys()
        if 'add_ACE' in keys:
            target_path = edit['add_ACE'].get('name', 'ACE')
            from_path = edit['add_ACE'].get('from')
            self.add_ACE(new_model, target_path, from_path)
        elif 'add_NME' in keys:
            target_path = edit['add_NME'].get('name', 'NME')
            from_path = edit['add_NME'].get('from')
            self.add_NME(new_model, target_path, from_path)

        self.output_model = new_model

        return self

    def add_ACE(self, model, target_path, from_path):
        logger.info("add_ACE from: {}".format(from_path))

        base_res = self._select_atomgroup(self.model, from_path)
        modeling = bridge.Modeling()
        atomgroup_ace = modeling.get_ACE(base_res)

        atomgroup_ace.path = target_path

        self.insert_atomgroup(model, target_path, atomgroup_ace)
        self.output_model = model

    def add_NME(self, model, target_path, from_path):
        logger.info("add_NME from: {}".format(from_path))

        base_res = self._select_atomgroup(self.model, from_path)
        # print('>>>> base_res')
        # print(base_res)
        modeling = bridge.Modeling()
        atomgroup_nme = modeling.get_NME(base_res)
        # atomgroup_nme = bridge.AtomGroup()

        atomgroup_nme.path = target_path

        self.insert_atomgroup(model, target_path, atomgroup_nme)
        self.output_model = model

    def insert_atomgroup(self, model, target_path, subgroup):
        '''Set the subgroup to the target path for the input model
        '''
        assert(isinstance(model, bridge.AtomGroup))
        assert(isinstance(subgroup, bridge.AtomGroup))

        paths = target_path.split('/')
        paths = [x for x in paths if (len(x) > 0)]
        current = model
        for path in paths:
            if not current.has_groupkey(path):
                current.set_group(path, bridge.AtomGroup())
                assert(current.has_group(path))
            current = current.get_group(path)

        current.name = subgroup.name
        for atom in subgroup.get_atom_list():
            current.set_atom(atom.name, atom)
