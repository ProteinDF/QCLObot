#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import copy
import pprint

try:
    import msgpack
except:
    import msgpack_pure as msgpack

import logging
logger = logging.getLogger(__name__)

import proteindf_bridge as bridge
from .utils import check_format_model

class TaskObject(object):
    """Task management base class

    - create task directory
    name: Name of the task. It is used for the work directory.


    Example:
    >>> task1 = TaskObject(name = 'test_task1')
    >>> print(task1.name)
    test_task1
    >>> os.path.exists('task1')
    True

    >>> task2 = TaskObject(name = 'test_task2')
    >>> print(task2.name)
    test_task2
    >>> os.path.exists('task2')
    True

    >>> task1_copy = TaskObject(task1)
    >>> print(task1_copy.name)
    test_task1

    """

    def __init__(self, *args, **kwargs):
        """initialize object

        :param str name: task name
        """
        self._initialize()

        if len(args) == 1:
            if isinstance(args[0], TaskObject):
                self._copy_constructor(args[0])
            elif isinstance(args[0], str):
                self._data['name'] = str(args[0])
            else:
                logger.critical('type mismatch')
                raise

        if 'name' in kwargs.keys():
            self._data['name'] = str(kwargs.get('name'))

        # check
        if 'name' not in self._data:
            raise

        # directory
        self._prepare_work_dir()


    def __del__(self):
        self.save()


    def _initialize(self):
        self._data = {}
        self._data['state_filename'] = 'qclobot_state.mpac'

        # not stored data
        self._basedir = os.path.abspath(os.curdir)

    def _copy_constructor(self, rhs):
        assert(isinstance(rhs, TaskObject))
        self._data = copy.deepcopy(rhs._data)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------
    # name -------------------------------------------------------------
    def _get_name(self):
        return self._data.get('name')
    name = property(_get_name)


    # work_dir ---------------------------------------------------------
    def _get_work_dir(self):
        """get working directory

        todo: mangling path
        """
        work_dir = os.path.join(self._basedir, self.name)
        return work_dir
    work_dir = property(_get_work_dir)


    # state_filename ---------------------------------------------------
    def _get_state_filename(self):
        return self._data.get('state_filename')
    state_filename = property(_get_state_filename)


    # model (input atomgroup) -----------------------------------------
    #   "model" data is an AtomGroup object formatted by 'MODEL', not 'protein'.
    def _get_model(self):
        answer = None
        model_raw_data = self._data.get("model", None)
        if model_raw_data != None:
            answer = bridge.AtomGroup(model_raw_data)
        return answer
    def _set_model(self, model):
        if check_format_model(model):
            self._data['model'] = model.get_raw_data()
        else:
            logger.critical("not support the format; use model format")
            raise
    model = property(_get_model, _set_model)


    # output_model -----------------------------------------------------
    def _get_output_model(self):
        answer = None
        model_raw_data = self._data.get("output_model", None)
        if model_raw_data != None:
            answer = bridge.AtomGroup(model_raw_data)
        return answer
    def _set_output_model(self, model):
        assert(check_format_model(model))
        self._data['output_model'] = model.get_raw_data()
    output_model = property(_get_output_model, _set_output_model)



    # ------------------------------------------------------------------
    # Archive
    # ------------------------------------------------------------------
    def load(self):
        path = os.path.join(self.work_dir, self.state_filename)
        if os.path.exists(path):
            logger.debug('load the fragment state: {}'.format(path))
            with open(path, 'rb') as f:
                packed = f.read()
                state_dat = msgpack.unpackb(packed)
                state_dat = bridge.Utils.to_unicode_dict(state_dat)
                self.set_by_raw_data(state_dat)
        else:
            logger.debug('not found the state file')

    def save(self):
        path = os.path.join(self.work_dir, self.state_filename)
        logger.debug('save the fragment state: {}'.format(path))

        state_dat = self.get_raw_data()
        # pprint.pprint(state_dat)
        packed = msgpack.packb(state_dat)
        with open(path, 'wb') as f:
            f.write(packed)

    def get_raw_data(self):
        return self.__getstate__()

    def set_by_raw_data(self, raw_data):
        self.__setstate__(raw_data)

    def __getstate__(self):
        state = copy.deepcopy(self._data)
        return state

    def __setstate__(self, state):
        assert(isinstance(state, dict))

        self._initialize()
        self._data['name'] = state.get('name')
        self._data['state_filename'] = state.get('state_filename')

    # ------------------------------------------------------------------
    # Work dir
    # ------------------------------------------------------------------
    def _prepare_work_dir(self):
        '''
        make working directory which is called its "name" property.
        '''
        assert(len(self.name) > 0)

        if not os.path.exists(self.work_dir):
            logger.debug('make work dir: {}'.format(self.work_dir))
            os.mkdir(self.work_dir)
            self.save()
        else:
            logger.debug('make work dir, but it already exists: {}'.format(self.work_dir))
            self.load()


    def cd_workdir(self, job_name=''):
        '''
        作業ディレクトリをオブジェクトのwork_dirに移動する
        '''
        logger.info('=' * 20)
        logger.info('>>>> {job_name}@{frame_name}'.format(job_name = job_name,
                                                          frame_name = self.name))
        logger.info('work dir: {work_dir}'.format(work_dir=self.work_dir))

        logger.info('=' * 20)
        os.chdir(self.work_dir)


    def restore_cwd(self):
        '''
        base ディレクトリに戻す
        '''
        os.chdir(self._basedir)
        logger.info('<<<< (basedir: {})\n'.format(self._basedir))


    def write_output_model(self, output_path):
        self.atomgroup2file(self.output_model, output_path)


    def atomgroup2file(self, atomgroup, output_path):
        '''output_pathの拡張子に応じて、bridge形式またはpdb形式でatomgroupをファイルに書き出す
        '''
        assert(isinstance(atomgroup, bridge.AtomGroup))

        abspath = os.path.abspath(output_path)
        (basename, ext) = os.path.splitext(abspath)
        ext = ext.lower()
        if ext in (".pdb", ".ent"):
            self.atomgroup2pdb(atomgroup, abspath)
        else:
            logger.info("save {path} as bridge file.".format(path=abspath))
            with open(abspath, "wb") as f:
                mpac_data = msgpack.packb(atomgroup.get_raw_data())
                f.write(mpac_data)


    def atomgroup2pdb(self, atomgroup, pdbfile,
                      model_name="model_1"):
        '''atomgroupをpdb形式で出力する

        atomgroupがmodels(複数のmodelで構成されている)の場合はそのまま出力する。
        atomgroupがmodelの場合は、model_nameを付加して出力する。
        '''
        assert(isinstance(pdbfile, str))

        self.cd_workdir()
        pdb = bridge.Pdb(mode = 'amber')

        protein = atomgroup
        if check_format_model(atomgroup):
            # transform MODEL object to the protein(models)
            # which has only one model.
            protein = bridge.AtomGroup()
            protein.set_group(model_name, atomgroup)

        pdb.set_by_atomgroup(protein)
        with open(pdbfile, 'w') as f:
            f.write(str(pdb))

        self.restore_cwd()


if __name__ == '__main__':
    #import sys,os
    #sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    import doctest
    doctest.testmod()
