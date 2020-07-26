#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import proteindf_bridge as bridge

from .modeler_task_object import ModelerTaskObject
from .process import Process

import logging
logger = logging.getLogger(__name__)


class QcProtonate(ModelerTaskObject):
    ''' execute protonate

    >>> tmp_pdb = bridge.Pdb('./data/sample/1AKG.pdb')
    >>> atomgroup = tmp_pdb.get_atomgroup()
    >>> p = QcProtonate(name='protonate_1akg', atomgroup=atomgroup)
    >>> p.protonate()
    USER  MOD reduce.3.24.130724 H: found=0, std=0, add=98, rem=0, adj=4
    ...

    >>> p.protonate_group()
    '''

    def __init__(self, parent, task):
        """ initialize protonate object

        :param str pdbfile: pdb file for protonation
        """
        super().__init__(parent, task)

        self._AMBERHOME = os.environ.get('AMBERHOME', '')
        self._check_AMBERHOME()

    def _check_AMBERHOME(self):
        if len(self._AMBERHOME) == 0:
            logger.warning("environ parameter, AMBERHOME, looks like empty.")

    ####################################################################
    # property
    ####################################################################

    # backend ----------------------------------------------------------
    def _get_backend(self):
        if 'backend' not in self._data:
            self._data['backend'] = 'reduce'
        return self._data.get('backend')
    backend = property(_get_backend)

    # model_name -------------------------------------------------------
    def _get_model_name(self):
        return self._data.get("model_name", "model_1")
    model_name = property(_get_model_name)

    ####################################################################
    # method
    ####################################################################

    def run(self, output_path=""):
        return_code = -1
        self.cd_workdir()

        if self.backend == 'reduce':
            input_pdbfile = os.path.join(self.work_dir, 'original.pdb')
            self.atomgroup2pdb(self.model, input_pdbfile,
                               model_name=self.model_name)

            out_pdbfile = os.path.join(self.work_dir, 'protonated.pdb')
            return_code = self._run_reduce(input_pdbfile,
                                           out_pdbfile)
        if return_code == 0:
            output_atomgroup = self._pdb2brd(out_pdbfile)

            # pickup first model as result
            self.output_model = output_atomgroup.get_group(self.model_name)

            # if len(output_path) > 0:
            #     output_path = os.path.join(self.work_dir, output_path)
            #     logger.info("output protonated file: {}".format(output_path))
            #     protein = bridge.AtomGroup()
            #     protein.set_group("model_1", self.output_model)
            #     self.atomgroup2file(protein, output_path)

        self.restore_cwd()
        return return_code

    def _pdb2brd(self, pdbfile):
        assert(isinstance(pdbfile, str))

        logger.info('pdb2brd: from {}'.format(pdbfile))
        pdb = bridge.Pdb(pdbfile, mode='amber')
        return pdb.get_atomgroup()

    def _run_reduce(self, in_pdbfile, out_pdbfile):
        assert(isinstance(in_pdbfile, str))
        assert(isinstance(out_pdbfile, str))

        p = Process()
        reduce_cmd = os.path.join(self._AMBERHOME, 'bin', 'reduce')
        cmd = "{} {}".format(reduce_cmd,
                             in_pdbfile)
        p.cmd(cmd)
        return_code = p.commit(out_pdbfile,
                               stdout_through=False,
                               stderr_through=False)

        return return_code

    def protonate_group(self):
        d_atomgroup = self.output_model ^ self.model

        d_path = os.path.join(self.work_dir, 'add_group.brd')
        bridge.save_msgpack(d_atomgroup.get_raw_data(), d_path)

    ####################################################################
    # Archive
    ####################################################################
    def __setstate__(self, state):
        super(QcProtonate, self).__setstate__(state)

        if "backend" in state:
            self._data["backend"] = state["backend"]
        if "model_name" in state:
            self._data["model_name"] = state["model_name"]


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
