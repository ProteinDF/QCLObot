#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from .taskobject import TaskObject
import proteindf_bridge as bridge


class MdObject(object):
    """Operate MD task

    : Example
    >>> md1 = MdObject(model)
    """

    def __init__(self, model, work_dir):
        assert(isinstance(model, bridge.AtomGroup))
        self._initialize()

        self.model = model
        self.work_dir = work_dir

    def _initialize(self):
        self._data = {}
        self._output_model = None

    # -------------------------------------------------------------------------
    # properties
    # -------------------------------------------------------------------------

    # work_dir
    def _get_work_dir(self):
        return self._work_dir

    def _set_work_dir(self, path):
        self._work_dir = str(path)

    work_dir = property(_get_work_dir, _set_work_dir)

    # model
    def _get_model(self):
        return self._model

    def _set_model(self, model):
        self._model = bridge.AtomGroup(model)
    model = property(_get_model, _set_model)

    # output_model
    def _get_output_model(self):
        return self._output_model

    def _set_output_model(self, model):
        assert(isinstance(model, bridge.AtomGroup))
        self._output_model = model
    output_model = property(_get_output_model, _set_output_model)

    # -------------------------------------------------------------------------
    # operation
    # -------------------------------------------------------------------------

    def opt(self):
        pass

    def md(self):
        pass


if __name__ == '__main__':
    import doctest
    doctest.testmod()
