#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import logging
logger = logging.getLogger(__name__)

try:
    import msgpack
except:
    import msgpack_pure as msgpack

import pdfbridge
from .taskobject import TaskObject
from .process import Process

class QcProtonate(TaskObject):
    ''' execute protonate

    >>> tmp_pdb = pdfbridge.Pdb('./data/sample/1AKG.pdb')
    >>> atomgroup = tmp_pdb.get_atomgroup()
    >>> p = QcProtonate(name='protonate_1akg', atomgroup=atomgroup)
    >>> p.protonate()
    USER  MOD reduce.3.24.130724 H: found=0, std=0, add=98, rem=0, adj=4
    ...

    >>> p.protonate_group()
    '''
    def __init__(self, name, atomgroup, backend='reduce'):
        """ initialize protonate object

        :param str pdbfile: pdb file for protonation
        """
        # initialize base object
        super().__init__(name=name)

        # initialize
        self._data['input_atomgroup'] = pdfbridge.AtomGroup(atomgroup)

        # backend
        self._data['backend'] = str(backend)
        self._AMBERHOME = os.environ.get('AMBERHOME', '')
        self._check_AMBERHOME()

        
    def _check_AMBERHOME(self):
        if len(self._AMBERHOME) == 0:
            logger.warning("environ parameter, AMBERHOME, looks like empty.")
        

    ####################################################################
    # property
    ####################################################################
    # input_atomgroup --------------------------------------------------
    def _get_input_atomgroup(self):
        return self._data.get('input_atomgroup')
    input_atomgroup = property(_get_input_atomgroup)

    # output_atomgroup -------------------------------------------------
    def _get_output_atomgroup(self):
        return self._data.get('output_atomgroup')
    output_atomgroup = property(_get_output_atomgroup)

    # backend ----------------------------------------------------------
    def _get_backend(self):
        return self._data.get('backend')
    backend = property(_get_backend)

    ####################################################################
    # method
    ####################################################################
    # run --------------------------------------------------------------
    def protonate(self):
        return_code = -1

        if self.backend == 'reduce':
            input_pdbfile = os.path.join(self.work_dir, 'original.pdb')
            self._brd2pdb(self.input_atomgroup, input_pdbfile)

            out_pdbfile=os.path.join(self.work_dir, 'protonated.pdb')
            return_code = self._run_reduce(input_pdbfile,
                                           out_pdbfile)
            if return_code == 0:
                output_atomgroup = self._pdb2brd(out_pdbfile)
                self._data['output_atomgroup'] = output_atomgroup
                output_brdfile = os.path.join(self.work_dir, 'protonated.brd')
                with open(output_brdfile, 'wb') as f:
                    raw_data = output_atomgroup.get_raw_data()
                    f.write(msgpack.packb(raw_data))

        return return_code

    def _brd2pdb(self, atomgroup, out_pdbfile):
        assert(isinstance(out_pdbfile, str))

        logger.info('brd2pdb: to {}'.format(out_pdbfile))
        pdb = pdfbridge.Pdb(mode = 'amber')
        pdb.set_by_atomgroup(atomgroup)
        with open(out_pdbfile, 'w') as f:
            f.write(str(pdb))

    def _pdb2brd(self, pdbfile):
        assert(isinstance(pdbfile, str))

        logger.info('pdb2brd: from {}'.format(pdbfile))
        pdb = pdfbridge.Pdb(pdbfile, mode = 'amber')
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
        d_atomgroup = self.output_atomgroup ^ self.input_atomgroup

        d_path = os.path.join(self.work_dir, 'add_group.brd')
        with open(d_path, 'wb') as f:
            raw_data = d_atomgroup.get_raw_data()
            f.write(msgpack.packb(raw_data))


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
