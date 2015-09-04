#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path
import tempfile
import logging

import qclobot as qclo

class QcProtonate(object):
    ''' execute protonate
    >>> p = QcProtonate()
    <BLANKLINE>
    >>> p.protonate('../data/sample/1AKG.pdb', '1AKG.addH.pdb')
    '''
    def __init__(self, pdbfile = ''):
        self._logger = logging.getLogger(__name__)
        
        self._AMBERHOME = os.environ.get('AMBERHOME', '')
        self._pdbfile = str(pdbfile)

        self._workdir = ''
        self._create_workdir()
        
    def __dell__(self):
        self._remove_workdir()
        
    # setup ------------------------------------------------------------
    def _create_workdir(self):
        self._workdir = tempfile.mkdtemp()

    def _remove_workdir(self):
        os.removedirs(self._WORKDIR)

    ####################################################################
    # property
    ####################################################################
    
    # pdbfile ----------------------------------------------------------
    def _set_pdbfile(self, pdbfile):
        self._pdbfile = str(pdbfile)

    def _get_pdbfile(self):
        return self._pdbfile

    pdbfile = property(_get_pdbfile,
                       _set_pdbfile)

    # run --------------------------------------------------------------
    def protonate(self, input_file, output_file):
        p = qclo.Process()

        reduce_cmd = os.path.join(self._AMBERHOME, 'bin', 'reduce')
        cmd = "{} {} {}".format(reduce_cmd,
                                input_file,
                                output_file)
        p.cmd(cmd)
        return_code = p.commit()

        
    def _run_reduce(self, output_file):
        p = qclo.Process()

        reduce_cmd = os.path.join(self._AMBERHOME, 'bin', 'reduce')
        cmd = "{} {} {}".format(reduce_cmd,
                          self.pdbfile,
                          output_file)
        p.cmd(cmd)
        return_code = p.commit()
        
    

if __name__ == '__main__':
    import doctest
    doctest.testmod()
