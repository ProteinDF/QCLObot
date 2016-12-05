#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import copy

try:
    import msgpack
except:
    import msgpack_pure as msgpack

import logging
logger = logging.getLogger(__name__)

import pdfbridge


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

    def _initialize(self):
        self._data = {}
        self._data['state_filename'] = 'qclobot_state.mpac'

        # not stored data
        self._basedir = os.path.abspath(os.curdir)
        # self._prev_dir = ''
        
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
                state_dat = pdfbridge.Utils.to_unicode_dict(state_dat)
                self.set_by_raw_data(state_dat)
        else:
            logger.debug('not found the state file')

    def save(self):
        path = os.path.join(self.work_dir, self.state_filename)
        logger.debug('save the fragment state: {}'.format(path))

        state_dat = self.get_raw_data()
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
        # self._prev_dir = os.path.abspath(os.curdir)
        os.chdir(self.work_dir)

    def restore_cwd(self):
        '''
        base ディレクトリに戻す
        '''
        os.chdir(self._basedir)
        logger.info('<<<< (basedir: {})\n'.format(self._basedir))
        
    
if __name__ == '__main__':
    #import sys,os
    #sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    import doctest
    doctest.testmod()
