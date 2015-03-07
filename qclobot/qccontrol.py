#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2002-2014 The ProteinDF project
# see also AUTHORS and README.
#
# This file is part of ProteinDF.
#
# ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import yaml
import logging
import pprint
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import pdfbridge as bridge
import qclobot as qclo

class QcControl(object):
    _modeling = bridge.Modeling()
    
    def __init__(self):
        self._logger = logging.getLogger(__name__)
        self._senarios = []
        self._cache = {}
        
        self._frames = {}
        self._frames['default'] = {}
        self._frames['default']['basis_set'] = 'DZVP2'
        self._frames['default']['brd_file'] = ''
        
    def run(self, path):
        self._load_yaml(path)
        for senario in self._senarios:
            for frame in senario:
                self._exec_frame(frame)
        
    def _load_yaml(self, path):
        f = open(path)
        contents = f.read()
        f.close()
        contents = contents.decode('utf8')

        self._senarios = []
        for d in yaml.load_all(contents):
            self._senarios.append(d)
        # self._logger.debug(pprint.pformat(self._senarios))

    def _save_yaml(self, data, path):
        assert(isinstance(data, dict))
        f = open(path, 'w')
        yaml.dump(data, f, encoding='utf8', allow_unicode=True)
        f.close()


    # ------------------------------------------------------------------
    def _exec_frame(self, frame_data):
        assert(isinstance(frame_data, dict))

        # --------------------------------------------------------------
        # setup
        # --------------------------------------------------------------
        # name
        if 'name' not in frame_data:
            self._logger.critical('name is not defined.')
        frame_name = frame_data.get('name')        
        frame = qclo.QcFrame(frame_name)

        self._logger.info('--- FRAME: {} ---'.format(frame_name))
        
        # frame configuration
        frame_basis_set = frame_data.get('basis_set',
                                         self._frames['default']['basis_set'])
        if 'XC_functional' in frame_data:
            frame.XC_functional = frame_data.get('XC_functional')
        if 'J_engine' in frame_data:
            frame.J_engine = frame_data.get('J_engine')
        if 'K_engine' in frame_data:
            frame.K_engine = frame_data.get('K_engine')
        if 'XC_engine' in frame_data:
            frame.XC_engine = frame_data.get('XC_engine')
        
        # coord
        frame_brd_file = frame_data.get('brd_file',
                                        self._frames['default']['brd_file'])
            

        #fragments_data = frame_data.get('fragments', [])
        #
        #for frg_data in fragments_data:
        #    # default parameter
        #    if 'brd_file' not in frg_data:
        #        frg_data['brd_file'] = frame_brd_file
        #    if 'basis_set' not in frg_data:
        #        frg_data['basis_set'] = frame_basis_set
        #    
        #    fragment = None
        #    if 'add_ACE' in frg_data:
        #        fragment = self._frg_add_ACE(frg_data)
        #    elif 'add_NME' in frg_data:
        #        fragment = self._frg_add_NME(frg_data)
        #    else:
        #        fragment = self._frg_default(frg_data)
        #    assert(fragment != None)
        #    frame[fragment.name] = fragment

        # --------------------------------------------------------------
        # default
        # --------------------------------------------------------------
        if frame_name == 'default':
            self._frames['default']['basis_set'] = frame_basis_set
            self._frames['default']['brd_file'] = frame_brd_file
            return
            
        # fragments
        self._logger.info('::make fragments')
        self._set_default(self._frames['default'], frame_data)
        fragments_list = self._get_fragments(frame_data.get('fragments', []), frame_data)
        for f in fragments_list:
            assert(isinstance(f, qclo.QcFragment))
            frame[f.name] = f
            self._logger.info('  > fragment: {}, parent={}'.format(f.name, frame[f.name].qc_parent.name))
            for subgrp_name, subgrp in frame[f.name].groups():
                self._logger.info('  > subfragment: {}, parent={}'.format(subgrp_name, subgrp.qc_parent.name))
            
        # --------------------------------------------------------------
        # action
        # --------------------------------------------------------------
        self._logger.info('::action')
        self._frames[frame_name] = frame

        # pre-SCF
        if frame_data.get('pre_scf', False):
            frame.calc_preSCF()
        
        # guess
        self._logger.info('::GUESS')
        guess = frame_data.get('guess', 'harris')
        guess = guess.lower()
        if guess == 'density':
            frame.guess_density()
        elif guess == 'qclo':
            frame.guess_QCLO()
        else:
            pass
            #frame.guess_Harris()

        # Single point calc.
        if frame_data.get('sp', False):
            frame.calc_sp()

    # ------------------------------------------------------------------
    # fragment
    # ------------------------------------------------------------------
    def _get_fragments(self, fragments_data, default):
        assert(isinstance(fragments_data, list))

        fragment = qclo.QcFragment()
        answer = []
        for frg_data in fragments_data:
            assert(isinstance(frg_data, dict))

            self._set_default(default, frg_data)
            name = frg_data.get('name', '')

            subfrg = None
            if 'fragments' in frg_data:
                subfrg_list = self._get_fragments(frg_data.get('fragments'), default)
                subfrg = qclo.QcFragment(name=name)
                for item in subfrg_list:
                    subfrg[item.name] = item
            elif 'add_ACE' in frg_data:
                subfrg = self._get_add_ACE(frg_data)
            elif 'add_NME' in frg_data:
                subfrg = self._get_add_NME(frg_data)
            elif 'brd_select' in frg_data:
                subfrg = self._get_default_fragment(frg_data)
            elif 'reference' in frg_data:
                subfrg = self._get_reference_fragment(frg_data)
            else:
                raise
            answer.append(subfrg)

        return answer

    def _set_default(self, default_values, update_values):
        keywords = ['brd_file', 'basis_set']
        for keyword in keywords:
                update_values.setdefault(keyword, default_values.get(keyword, None))
        
    def _get_default_fragment(self, frg_data):
        assert(isinstance(frg_data, dict))
        
        atomgroup = None
        brd_select = frg_data.get('brd_select')
        brd_file_path = frg_data.get('brd_file')
        atomgroup = self._select_atomgroup(brd_file_path, brd_select)

        frg = qclo.QcFragment(atomgroup)
        frg.margin = False
        if 'name' not in frg_data:
            raise
        frg.name = frg_data.get('name')

        self._set_basis_set(frg, frg_data.get('basis_set'))
        
        return frg

    def _get_reference_fragment(self, frg_data):
        assert(isinstance(frg_data, dict))

        ref_frame = frg_data['reference']['frame']
        ref_fragment = frg_data['reference']['fragment']

        return self._frames[ref_frame][ref_fragment]
        
    def _get_add_ACE(self, frg_data):
        assert(isinstance(frg_data, dict))
        brd_file_path = frg_data.get('brd_file')
        brd_select = frg_data.get('brd_select')
        atomgroup = self._select_atomgroup(brd_file_path, brd_select)
        ag_ACE = self._modeling.get_ACE_simple(atomgroup)

        ACE = qclo.QcFragment(ag_ACE)
        ACE.margin = True
        ACE.name = 'ACE'

        self._set_basis_set(ACE, frg_data.get('basis_set'))

        return ACE
        
    def _get_add_NME(self, frg_data):
        assert(isinstance(frg_data, dict))
        brd_file_path = frg_data.get('brd_file')
        brd_select = frg_data.get('brd_select')
        atomgroup = self._select_atomgroup(brd_file_path, brd_select)
        ag_NME = self._modeling.get_NME_simple(atomgroup)
        
        NME = qclo.QcFragment(ag_NME)
        NME.margin = True
        NME.name = 'NME'

        self._set_basis_set(NME, frg_data.get('basis_set'))

        return NME

    def _select_atomgroup(self, brd_file_path, brd_select):
        assert(isinstance(brd_file_path, str))
        assert(isinstance(brd_select, str))
        # self._logger.debug('_get_atomgroup(): brd_path={}, select={}'.format(brd_file_path, brd_select))

        self._cache.setdefault('brdfile', {})
        if brd_file_path not in self._cache['brdfile']:
            brd_fh = open(brd_file_path, 'rb')
            brd_data = msgpack.unpackb(brd_fh.read())
            brd_fh.close()

            self._cache['brdfile'][brd_file_path] = {}
            self._cache['brdfile'][brd_file_path]['atomgroup'] = bridge.AtomGroup(brd_data)
        
        atomgroup = self._cache['brdfile'][brd_file_path]['atomgroup']
        selecter = bridge.Select_Path(brd_select)
        answer = atomgroup.select(selecter)

        return answer

    def _set_basis_set(self, obj, basis_set_name):
        assert(isinstance(basis_set_name, str))

        if isinstance(obj, qclo.QcFragment):
            for subfrg_name, subfrg in obj.groups():
                self._set_basis_set(subfrg, basis_set_name)
            for atm_name, atm in obj.atoms():
                atm.basisset = '{}.{}'.format(basis_set_name, atm.symbol)
        else:
            raise
        
        
if __name__ == '__main__':
    pass
    
