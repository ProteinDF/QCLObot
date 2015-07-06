#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2002-2015 The ProteinDF project
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
import jinja2
    
import pdfbridge as bridge
import qclobot as qclo

class QcControl(object):
    _modeling = bridge.Modeling()
    
    def __init__(self):
        self._logger = logging.getLogger(__name__)
        self._senarios = []
        self._cache = {}

        self._vars = {}
        
        self._frames = {}
        self._frames['default'] = {}
        self._frames['default']['basis_set'] = 'DZVP2'
        self._frames['default']['brd_file'] = ''
        
    def run(self, path):
        self._load_yaml(path)

        # exec senarios
        for senario in self._senarios:
            if 'vars' in senario:
                self._run_vars(senario['vars'])
            if 'tasks' in senario:
                tasks = senario['tasks']
                for task in tasks:
                    self._run_frame(task)
        
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
    # vars
    # ------------------------------------------------------------------
    def _run_vars(self, in_vars_data):
        assert(isinstance(in_vars_data, dict))

        self._vars = dict(in_vars_data)
        
    # ------------------------------------------------------------------
    # task or frame
    # ------------------------------------------------------------------
    def _run_frame(self, frame_data):
        assert(isinstance(frame_data, dict))
        
        # condition ----------------------------------------------------
        if self._exec_condition_with_items(frame_data):
            return
        if self._exec_condition_when(frame_data):
            return
        
        # task ---------------------------------------------------------
        self._exec_task_debug(frame_data)

        if self._exec_task_mail(frame_data):
            return

        # frame --------------------------------------------------------
        if self._exec_frame_default(frame_data):
            return

        self._exec_frame_object(frame_data)
        
    # ------------------------------------------------------------------
    # condition
    # ------------------------------------------------------------------
    def _exec_condition_with_items(self, frame_data):
        is_break = False

        with_items = frame_data.get('with_items', None)
        if with_items:
            iter_items = list(self._vars.get(str(with_items), []))

            new_frame_data = dict(frame_data)
            new_frame_data.pop('with_items')
            yaml_str = yaml.dump(new_frame_data)
            template = jinja2.Template(yaml_str)

            for item in iter_items:
                yaml_str = template.render(item = item)
                new_frame_data = yaml.load(yaml_str)
                self._run_frame(new_frame_data)
            is_break = True

        return is_break

    def _exec_condition_when(self, task_data):
        is_break = False

        when_phrase = task_data.get('when', None)
        if when_phrase:
            globals_data = self._vars
            judge = eval(when_phrase, globals_data)
            if judge:
                new_task_data = dict(task_data)
                new_task_data.pop('when')
                self._run_frame(new_task_data)
            is_break = True
        
        return is_break

    # ------------------------------------------------------------------
    # task
    # ------------------------------------------------------------------
    def _exec_task_debug(self, in_task_data):
        is_break = False

        return is_break
        
    def _exec_task_mail(self, in_task_data):
        is_break = False

        if 'mail' in in_task_data:
            mail_data = self._get_value('mail', in_task_data)
            mailer = bridge.Mail()
            mailer.smtp_server = mail_data.get('smtp_server')
            if 'smtp_port' in mail_data:
                mailer.smtp_port = mail_data.get('smtp_port')
            if 'use_SSL' in mail_data:
                mailer.use_SSL = mail_data.get('use_SSL')
            mailer.smtp_account = mail_data.get('smtp_account')
            mailer.smtp_password = mail_data.get('smtp_password')
            mailer.from_address = mail_data.get('from_address')
            mailer.to_address = mail_data.get('to_address')
            mailer.subject = mail_data.get('subject')
            mailer.text = mail_data.get('msg')
            mailer.send()

            is_break = True
            
        return is_break
        
    # ------------------------------------------------------------------
    # frame
    # ------------------------------------------------------------------
    def _exec_frame_default(self, in_frame_data):
        assert(isinstance(in_frame_data, dict))
        is_break = False
        
        name = in_frame_data.get('name', None)
        if ((name != None) and (name == 'default')):
            default_data = dict(in_frame_data)
            default_data.pop('name')
            self._frames['default'].update(default_data)
            is_break = True

        return is_break

    def _exec_frame_object(self, frame_data):
        assert(isinstance(frame_data, dict))
        
        # --------------------------------------------------------------
        # setup
        # --------------------------------------------------------------
        # name
        if 'name' not in frame_data:
            self._logger.critical('name is not defined.')
        frame_name = frame_data.get('name')        
        self._logger.info('--- FRAME: {} ---'.format(frame_name))
        frame = qclo.QcFrame(frame_name)

        # setup default value
        self._set_default(self._frames['default'], frame_data)

        # frame configuration
        charge = self._get_value('charge', frame_data)
        if charge:
            frame.charge = charge
        XC_functional = self._get_value('XC_functional', frame_data)
        if XC_functional:
            frame.pdfparam.xc_functional = XC_functional
        J_engine = self._get_value('J_engine', frame_data)
        if J_engine:
            frame.pdfparam.j_engine = J_engine
        K_engine = self._get_value('K_engine', frame_data)
        if K_engine:
            frame.pdfparam.k_engine = K_engine
        XC_engine = self._get_value('XC_engine', frame_data)
        if XC_engine:
            frame.pdfparam.xc_engine = XC_engine

        # fragments
        self._logger.info('::make fragments')
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
    # utils
    # ------------------------------------------------------------------
    def _get_value(self, keyword, input_data):
        '''
        return a value corresponding to the keyword from the input_data or defaults
        '''
        assert(isinstance(keyword, str))

        answer = None
        if keyword in input_data:
            answer = input_data.get(keyword)
        elif keyword in self._frames['default']:
            answer = self._frames['default'][keyword]
        return answer
        
    # ------------------------------------------------------------------
    # fragment
    # ------------------------------------------------------------------
    def _get_fragments(self, fragments_data, default):
        assert(isinstance(fragments_data, list))
        self._logger.debug(">>>> _get_fragments()")
        self._logger.debug(str(default))
        
        fragment = qclo.QcFragment()
        answer = []
        for frg_data in fragments_data:
            assert(isinstance(frg_data, dict))

            self._set_default(default, frg_data)
            self._logger.debug(">>>> _get_fragments(): loop")
            self._logger.debug(frg_data)
            name = frg_data.get('name', '')

            subfrg = None
            if 'fragments' in frg_data:
                subfrg_list = self._get_fragments(frg_data.get('fragments'), default)
                subfrg = qclo.QcFragment(name=name)
                for item in subfrg_list:
                    if item == None:
                        self._logger.debug('name: {}'.format(name))
                        self._logger.debug('subfrg_list: {}'.format(str(subfrg_list)))
                        raise qclo.QcControlError('unknown subfragments:', str(frg_data))
                    subfrg[item.name] = item
            elif 'add_H' in frg_data:
                subfrg = self._get_add_H(frg_data)
            elif 'add_CH3' in frg_data:
                subfrg = self._get_add_CH3(frg_data)
            elif 'add_ACE' in frg_data:
                subfrg = self._get_add_ACE(frg_data)
            elif 'add_NME' in frg_data:
                subfrg = self._get_add_NME(frg_data)
            elif 'brd_select' in frg_data:
                subfrg = self._get_default_fragment(frg_data)
            elif 'reference' in frg_data:
                subfrg = self._get_reference_fragment(frg_data)
            else:
                self.critical(str(frg_data))
                raise qclo.QcControlError('unknown fragment:', str(frg_data))

            answer.append(subfrg)

        return answer

    def _set_default(self, default_values, update_values):
        '''
        copy default values to variable
        '''
        keywords = ['brd_file',
                    'basis_set',
                    'guess',
                    'XC_functional',
                    'J_engine',
                    'K_engine',
                    'XC_engine']
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

    def _get_add_H(self, frg_data):
        assert(isinstance(frg_data, dict))
        brd_file_path = frg_data.get('brd_file')
        brd_select_H = frg_data.get('displacement')
        atomgroup_H = self._select_atomgroup(brd_file_path, brd_select_H)
        atomgroup_H = atomgroup_H.get_atomlist()
        assert(atomgroup_H.get_number_of_atoms() > 0)
        (key_H, atom_H) = list(atomgroup_H.atoms())[0]
        atom_H.symbol = 'H'
        
        ag_H = bridge.AtomGroup()
        ag_H.set_atom('H*', atom_H)
        
        H = qclo.QcFragment(ag_H)
        H.margin = True
        if 'name' not in frg_data:
            raise
        H.name = frg_data.get('name')

        self._set_basis_set(H, frg_data.get('basis_set'))

        return H
        
    def _get_add_CH3(self, frg_data):
        assert(isinstance(frg_data, dict))
        brd_file_path = frg_data.get('brd_file')
        brd_select_C1 = frg_data.get('displacement')
        brd_select_C2 = frg_data.get('root')

        atomgroup_C1 = self._select_atomgroup(brd_file_path, brd_select_C1)
        atomgroup_C1 = atomgroup_C1.get_atomlist()
        assert(atomgroup_C1.get_number_of_atoms() > 0)
        (key_C1, atom_C1) = list(atomgroup_C1.atoms())[0]

        atomgroup_C2 = self._select_atomgroup(brd_file_path, brd_select_C2)
        atomgroup_C2 = atomgroup_C2.get_atomlist()
        assert(atomgroup_C2.get_number_of_atoms() > 0)
        (key_C2, atom_C2) = list(atomgroup_C2.atoms())[0]
        ag_CH3 = self._modeling.add_methyl(atom_C1, atom_C2)
        ag_CH3.set_atom('C', atom_C1)
        
        CH3 = qclo.QcFragment(ag_CH3)
        CH3.margin = True
        if 'name' not in frg_data:
            raise
        CH3.name = frg_data.get('name')

        self._set_basis_set(CH3, frg_data.get('basis_set'))

        return CH3
        
    def _get_add_ACE(self, frg_data):
        assert(isinstance(frg_data, dict))
        brd_file_path = frg_data.get('brd_file')
        brd_select = frg_data.get('brd_select')
        atomgroup = self._select_atomgroup(brd_file_path, brd_select)
        ag_ACE = self._modeling.get_ACE_simple(atomgroup)

        ACE = qclo.QcFragment(ag_ACE)
        ACE.margin = True
        if 'name' not in frg_data:
            raise
        ACE.name = frg_data.get('name')

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
        if 'name' not in frg_data:
            raise
        NME.name = frg_data.get('name')
        # NME.name = 'NME'

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
                atm.basisset = 'O-{}.{}'.format(basis_set_name, atm.symbol)
                atm.basisset_j = 'A-{}.{}'.format(basis_set_name, atm.symbol)
                atm.basisset_xc = 'A-{}.{}'.format(basis_set_name, atm.symbol)
                atm.basisset_gridfree = 'O-{}.{}'.format(basis_set_name, atm.symbol)
        else:
            raise
        
        
if __name__ == '__main__':
    pass
    
