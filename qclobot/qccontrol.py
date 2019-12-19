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

import jinja2
import pprint
import yaml
import logging
logger = logging.getLogger(__name__)
try:
    import msgpack
except:
    import msgpack_pure as msgpack

from .qcerror import QcControlError
from .qccontrolobject import QcControlObject
from .qcfragment import QcFragment
from .qcframe import QcFrame
from . import __version__
import proteindf_bridge as bridge


class QcControl(QcControlObject):
    _modeling = bridge.Modeling()

    def __init__(self):
        super(QcControl, self).__init__()
        self._cache = {}

        self._frames = {}
        self._frames['default'] = {}
        self._frames['default']['basis_set'] = 'DZVP2'
        self._frames['default']['brd_file'] = ''

        self._last_frame_name = ''

    # ------------------------------------------------------------------
    # property
    # ------------------------------------------------------------------
    @property
    def last_frame_name(self):
        ''' return the last frame name
        '''
        return self._last_frame_name

    def get_frame(self, name):
        return self._frames.get(name, None)

    # ------------------------------------------------------------------
    # task or frame
    # ------------------------------------------------------------------
    def _run_task_cmd(self, frame_data):
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
            iter_items = list(self._vars.get(
                bridge.Utils.to_unicode(with_items), []))

            new_frame_data = dict(frame_data)
            new_frame_data.pop('with_items')
            yaml_str = yaml.dump(new_frame_data)
            template = jinja2.Template(yaml_str)

            for item in iter_items:
                logger.info('template render: item={}'.format(repr(item)))
                yaml_str = template.render(item=item)
                for new_frame_data in yaml.load_all(yaml_str, Loader=yaml.SafeLoader):
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
        '''
        execute something to frame object
        '''

        assert(isinstance(frame_data, dict))

        # --------------------------------------------------------------
        # setup
        # --------------------------------------------------------------
        # name
        if 'name' not in frame_data:
            logger.critical('name is not defined.')
        frame_name = frame_data.get('name')
        logger.info("#" * 72)
        logger.info('# FRAME: {} '.format(frame_name))
        logger.info("#" * 72)
        frame = QcFrame(frame_name)

        # setup default value
        self._set_default(self._frames['default'], frame_data)
        # print(">>>> _exec_frame_object")
        # print(self._frames['default'])
        # print(frame_data)
        # print("<<<<")

        # frame configuration
        frame.pdfparam.guess = self._get_value('guess', frame_data)

        charge = self._get_value('charge', frame_data)
        if charge:
            frame.charge = charge

        frame.pdfparam.CDAM_tau = self._get_value('CDAM_tau', frame_data)
        frame.pdfparam.CD_epsilon = self._get_value('CD_epsilon', frame_data)

        frame.pdfparam.convergence_threshold_energy = self._get_value(
            'convergence/threshold_energy', frame_data)
        frame.pdfparam.convergence_threshold = self._get_value(
            'convergence/threshold', frame_data)
        frame.pdfparam.convergence_type = self._get_value(
            'convergence/type', frame_data)

        frame.pdfparam.orbital_independence_threshold = self._get_value(
            'orbital_independence_threshold', frame_data)
        frame.pdfparam.orbital_independence_threshold_canonical = self._get_value(
            'orbital_independence_threshold/canonical', frame_data)
        frame.pdfparam.orbital_independence_threshold_lowdin = self._get_value(
            'orbital_independence_threshold/lowdin', frame_data)

        frame.pdfparam.scf_acceleration = self._get_value(
            'scf_acceleration', frame_data)
        frame.pdfparam.scf_acceleration_damping_start_number = self._get_value(
            'scf_acceleration/damping/start_number', frame_data)
        frame.pdfparam.scf_acceleration_damping_damping_factor = self._get_value(
            'scf_acceleration/damping/damping_factor', frame_data)
        frame.pdfparam.scf_acceleration_anderson_start_number = self._get_value(
            'scf_acceleration/anderson/start_number', frame_data)
        frame.pdfparam.scf_acceleration_anderson_damping_factor = self._get_value(
            'scf_acceleration/anderson/damping_factor', frame_data)

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

        frame.pdfparam.gridfree_dedicated_basis = self._get_value(
            'gridfree/dedicated_basis', frame_data)
        frame.pdfparam.gridfree_orthogonalize_method = self._get_value(
            'gridfree/orthogonalize_method', frame_data)
        frame.pdfparam.gridfree_CDAM_tau = self._get_value(
            'gridfree/CDAM_tau', frame_data)
        frame.pdfparam.gridfree_CD_epsilon = self._get_value(
            'gridfree/CD_epsilon', frame_data)

        if self._get_value('pdf_extra_keywords', frame_data) != None:
            frame.pdfparam.extra_keywords = self._get_value(
                'pdf_extra_keywords', frame_data)
                
        #print(">>>> qccontrol:")
        # print(repr(frame.pdfparam.extra_keywords))
        # print("<<<<")

        # cmd alias
        if self._get_value('cmd_alias', frame_data) != None:
            frame.set_command_alias(self._get_value('cmd_alias', frame_data))

        # fragments
        logger.info('> make fragments for [{frame_name}]'.format(
            frame_name=frame_name))
        fragments_list = self._get_fragments(
            frame_data.get('fragments', []), frame_data)
        for f in fragments_list:
            assert(isinstance(f, QcFragment))
            frame[f.name] = f

            # fragment information
            #parent_frame_name = ""
            # if f.parent != None:
            #    parent_frame_name = f.parent.name
            # logger.info("> add fragment {ref_frame}/{ref_fragment} to {frame_name}".format(
            #    ref_frame=parent_frame_name,
            #    ref_fragment=f.name,
            #    frame_name=frame_name
            # ))
            # for subgrp_name, subgrp in frame[f.name].groups():
            #    logger.info('  > subfragment: {}, parent={}'.format(subgrp_name, subgrp.parent.name))
        logger.debug('> make fragments for [{frame_name}]: end'.format(
            frame_name=frame_name))

        # --------------------------------------------------------------
        # action
        # --------------------------------------------------------------
        logger.info('::action')
        frame.save()
        self._frames[frame_name] = frame
        self._last_frame_name = frame_name

        # pre-SCF
        if frame_data.get('pre_scf', False):
            frame.calc_preSCF()

        # make guess by QCLObot
        logger.info('::GUESS')
        guess = frame_data.get('guess', 'harris')
        guess = guess.lower()
        guess_force = frame_data.get("guess/force", "False")
        if isinstance(guess_force, str):
            guess_force = guess_force.lower()
            if (len(guess_force) > 0) and (guess_force in ("true", "yes", "1")):
                guess_force = True
            else:
                guess_force = False

        if guess == 'density':
            frame.guess_density(force=guess_force)
        elif guess == 'qclo':
            frame.guess_QCLO(force=guess_force)
        else:
            pass
            # frame.guess_Harris()

        # Single point calc.
        if frame_data.get('sp', False):
            frame.calc_sp()

        # force
        if frame_data.get('force', False):
            frame.calc_force()

        # summary
        summary_act = frame_data.get('summary', False)
        if isinstance(summary_act, dict):
            format_str = summary_act.get('format', None)
            filepath = summary_act.get('filepath', None)
            frame.summary(format_str=format_str,
                          filepath=filepath)
        elif isinstance(summary_act, bridge.basestring):
            frame.summary(format_str=summary_act)
        elif isinstance(summary_act, bool):
            if summary_act:
                frame.summary()

        # finalize
        frame.save()

    # ------------------------------------------------------------------
    # utils
    # ------------------------------------------------------------------
    def _get_value(self, keyword, input_data):
        '''
        return a value corresponding to the keyword from the input_data or defaults
        '''
        keyword = bridge.Utils.to_unicode(keyword)

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
        logger.info("> make list of fragments: start")
        logger.debug(str(default))

        #fragment = QcFragment()
        answer = []
        for frg_data in fragments_data:
            assert(isinstance(frg_data, dict))

            self._set_default(default, frg_data)
            logger.debug("> get_fragments(): loop")
            logger.debug(frg_data)
            name = frg_data.get('name', '')
            logger.debug('> name: {}'.format(name))

            subfrg = None
            if 'fragments' in frg_data:
                logger.info(
                    "> create subfragment for {name}: start".format(name=name))
                subfrg_list = self._get_fragments(
                    frg_data.get('fragments'), default)
                subfrg = QcFragment(name=name)
                for item in subfrg_list:
                    logger.debug('> list name: {}'.format(item.name))
                    subfrg[item.name] = item
                logger.info(
                    "> create subfragment for {name}: end".format(name=name))
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
            elif 'atomlist' in frg_data:
                subfrg = self._getfrg_atomlist(frg_data)
            elif 'reference' in frg_data:
                subfrg = self._get_reference_fragment(frg_data)
            else:
                logger.critical(
                    "unknown fragment command/object: {}".format(str(frg_data.keys())))
                raise QcControlError(
                    'fragment', 'unknown fragment: {}'.format(str(frg_data)))

            if subfrg == None:
                frg_data_str = pprint.pformat(frg_data)
                logger.critical('> unknown sub-fragment:')
                logger.critical('>  sub-fragment information:')
                logger.critical(frg_data_str)
                raise QcControlError('unknown subfragment:', frg_data_str)

            assert(isinstance(subfrg, QcFragment))
            subfrg_atomgroup = subfrg.get_AtomGroup()
            logger.info("subfrg append: {name} ({formula}: elec={elec})".format(
                name=name,
                formula=subfrg_atomgroup.get_formula(),
                elec=subfrg_atomgroup.charge))
            answer.append(subfrg)

        logger.info("> make list of fragments: end")
        return answer

    def _set_default(self, default_values, update_values):
        '''
        copy default values to variable
        '''
        keywords = ['brd_file',
                    'basis_set',
                    'basis_set_aux',
                    'basis_set_gridfree',
                    'guess',
                    'cut_value',
                    'CDAM_tau',
                    'CD_epsilon',
                    'orbital_independence_threshold',
                    'orbital_independence_threshold/canonical',
                    'orbital_independence_threshold/lowdin',
                    'gridfree/dedicated_basis',
                    'gridfree/orthogonalize_method',
                    'gridfree/CDAM_tau',
                    'gridfree/CD_epsilon',
                    'scf_acceleration',
                    'scf_acceleration/damping/damping_factor',
                    'scf_acceleration/anderson/start_number',
                    'scf_acceleration/anderson/damping_factor',
                    'convergence/threshold_energy',
                    'convergence/threshold',
                    'XC_functional',
                    'J_engine',
                    'K_engine',
                    'XC_engine',
                    'pdf_extra_keywords']
        for keyword in keywords:
            update_values.setdefault(
                keyword, default_values.get(keyword, None))

    def _get_default_fragment(self, frg_data):
        assert(isinstance(frg_data, dict))

        atomgroup = None
        brd_select = frg_data.get('brd_select')
        brd_file_path = frg_data.get('brd_file')
        atomgroup = self._select_atomgroup(brd_file_path, brd_select)

        frg = QcFragment(atomgroup)
        frg.margin = False
        if 'name' not in frg_data:
            raise QcControlError('name keyword NOT FOUND', str(frg_data))
        frg.name = frg_data.get('name')

        self._set_basis_set(frg,
                            frg_data.get('basis_set'),
                            frg_data.get('basis_set_aux', None),
                            frg_data.get('basis_set_gridfree', None))

        return frg

    def _get_reference_fragment(self, frg_data):
        """return new QcFragment object referenced from input(frg_data)
        """
        assert(isinstance(frg_data, dict))

        if not 'frame' in frg_data['reference']:
            raise QcControlError('NOT FOUND frame key in reference fragment',
                                 pprint.pformat(frg_data))
        if not 'fragment' in frg_data['reference']:
            raise QcControlError('NOT FOUND fragment key in reference fragment',
                                 pprint.pformat(frg_data))

        ref_frame_name = bridge.Utils.to_unicode(
            frg_data['reference']['frame'])
        ref_fragment_name = bridge.Utils.to_unicode(
            frg_data['reference']['fragment'])

        if ref_frame_name not in self._frames:
            raise QcControlError('UNKNOWN FRAME',
                                 ref_frame_name)
        ref_frame = self._frames.get(ref_frame_name)

        if not ref_frame.has_fragment(ref_fragment_name):
            logger.critical(str(self._frames[ref_frame_name]))
            raise QcControlError('UNKNOWN FRAGMENT',
                                 '{} in {}'.format(ref_fragment_name,
                                                   ref_frame_name))

        logger.info(
            "> reference fragment: [{}/{}]".format(ref_frame_name, ref_fragment_name))
        ref_fragment = ref_frame[ref_fragment_name]

        fragment = QcFragment(ref_fragment)
        fragment.name = frg_data.get('name')
        fragment.ref_fragment = ref_fragment
        return fragment

    def _get_add_H(self, frg_data):
        assert(isinstance(frg_data, dict))
        brd_file_path = frg_data.get('brd_file')
        brd_select_H = frg_data.get('displacement')
        atomgroup_H = self._select_atomgroup(brd_file_path, brd_select_H)
        atomgroup_H = atomgroup_H.get_atom_list()
        assert(len(atomgroup_H) > 0)
        atom_H = atomgroup_H[0]
        atom_H.symbol = 'H'

        ag_H = bridge.AtomGroup()
        ag_H.set_atom('H*', atom_H)

        H = QcFragment(ag_H)
        H.margin = True
        if 'name' not in frg_data:
            raise QcControlError('NOT FOUND name key in add_H', str(frg_data))
        H.name = frg_data.get('name')

        self._set_basis_set(H,
                            frg_data.get('basis_set'),
                            frg_data.get('basis_set_aux', None),
                            frg_data.get('basis_set_gridfree', None))

        return H

    def _get_add_CH3(self, frg_data):
        assert(isinstance(frg_data, dict))
        brd_file_path = frg_data.get('brd_file')
        brd_select_C1 = frg_data.get('displacement')
        brd_select_C2 = frg_data.get('root')

        atomgroup_C1 = self._select_atomgroup(brd_file_path, brd_select_C1)
        atomgroup_C1 = atomgroup_C1.get_atom_list()
        if len(atomgroup_C1) == 0:
            logger.critical(
                'cannnot find displacement atom: key={}'.format(brd_select_C1))
            raise QcControlError('No atoms found for "displacement" in add_CH3',
                                 brd_select_C1)
        atom_C1 = atomgroup_C1[0]

        atomgroup_C2 = self._select_atomgroup(brd_file_path, brd_select_C2)
        atomgroup_C2 = atomgroup_C2.get_atom_list()
        if len(atomgroup_C2) == 0:
            logger.critical(
                'cannnot find displacement atom: key={}'.format(brd_select_C2))
            raise QcControlError('No atoms found for "root" in add_CH3')
        atom_C2 = atomgroup_C2[0]
        ag_CH3 = self._modeling.add_methyl(atom_C1, atom_C2)
        ag_CH3.set_atom('C', atom_C1)

        CH3 = QcFragment(ag_CH3)
        CH3.margin = True
        if 'name' not in frg_data:
            raise QcControlError(
                'NOT FOUND name key in add_CH3', str(frg_data))
        CH3.name = frg_data.get('name')

        self._set_basis_set(CH3,
                            frg_data.get('basis_set'),
                            frg_data.get('basis_set_gridfree', None))

        return CH3

    def _get_add_ACE(self, frg_data):
        assert(isinstance(frg_data, dict))
        brd_file_path = frg_data.get('brd_file')
        brd_select = frg_data.get('brd_select')
        atomgroup = self._select_atomgroup(brd_file_path, brd_select)
        ag_ACE = self._modeling.get_ACE_simple(atomgroup)

        ACE = QcFragment(ag_ACE)
        ACE.margin = True
        if 'name' not in frg_data:
            raise QcControlError(
                'NOT FOUND name key in add_ACE', str(frg_data))
        ACE.name = frg_data.get('name')

        self._set_basis_set(ACE,
                            frg_data.get('basis_set'),
                            frg_data.get('basis_set_aux', None),
                            frg_data.get('basis_set_gridfree', None))

        return ACE

    def _get_add_NME(self, frg_data):
        assert(isinstance(frg_data, dict))
        brd_file_path = frg_data.get('brd_file')
        brd_select = frg_data.get('brd_select')
        atomgroup = self._select_atomgroup(brd_file_path, brd_select)
        ag_NME = self._modeling.get_NME_simple(atomgroup)

        NME = QcFragment(ag_NME)
        NME.margin = True
        if 'name' not in frg_data:
            raise QcControlError(
                'NOT FOUND name key in add_NME', str(frg_data))
        NME.name = frg_data.get('name')

        self._set_basis_set(NME,
                            frg_data.get('basis_set'),
                            frg_data.get('basis_set_aux', None),
                            frg_data.get('basis_set_gridfree', None))

        return NME

    def _getfrg_atomlist(self, frg_data):
        assert(isinstance(frg_data, dict))
        atomlist = frg_data.get('atomlist')

        atomgroup = bridge.AtomGroup()
        if 'name' not in frg_data:
            raise QcControlError(
                'NOT FOUND name key in atomlist', str(frg_data))
        atomgroup.name = frg_data.get('name')

        index = 1
        if isinstance(atomlist, list):
            for line in atomlist:
                if isinstance(line, bridge.basestring):
                    if line[0] == '/':
                        # case '/model/...'
                        brd_select = line
                        brd_file_path = frg_data.get('brd_file')
                        selected_atomgroup = self._select_atomgroup(
                            brd_file_path, brd_select)
                        #logger.debug('select path: {}'.format(brd_select))
                        #logger.debug('selected: {}'.format(str(selected_atomgroup)))
                        if selected_atomgroup.get_number_of_all_atoms() > 0:
                            for selected_atom in selected_atomgroup.get_atom_list():
                                atomgroup.set_atom(index, selected_atom)
                                index += 1
                        else:
                            logger.warn("not math: {}".format(brd_select))
                        continue
                    else:
                        # case: "H 1.00 2.00 3.00"
                        line = line.split()
                if isinstance(line, list):
                    if len(line) == 4:
                        # case [H, 1.00, 2.00, 3.00]
                        atom = bridge.Atom()
                        atom.symbol = line[0]
                        atom.xyz = bridge.Position(line[1], line[2], line[3])
                        atomgroup.set_atom(index, atom)
                        index += 1
                    elif len(line) == 5:
                        # case [X, 1,00, 2.00, 3.00, -1.0] for dummy charge
                        atom = bridge.Atom()
                        atom.symbol = line[0]
                        atom.xyz = bridge.Position(line[1], line[2], line[3])
                        atom.charge = line[4]
                        atomgroup.set_atom(index, atom)
                        index += 1
                    else:
                        logger.error(
                            "mismatch atom data size: {}".format(str(line)))
                else:
                    logger.error(
                        "atomlist item shuld be array or string: {}".format(str(line)))
        else:
            logger.error("atomlist shuld be list: {}".format(str(atomlist)))

        frg = QcFragment(atomgroup)
        self._set_basis_set(frg,
                            frg_data.get('basis_set'),
                            frg_data.get('basis_set_aux', None),
                            frg_data.get('basis_set_gridfree', None))
        return frg

    def _select_atomgroup(self, brd_file_path, brd_select):
        brd_file_path = bridge.Utils.to_unicode(brd_file_path)
        brd_select = bridge.Utils.to_unicode(brd_select)
        # logger.debug('_get_atomgroup(): brd_path={}, select={}'.format(brd_file_path, brd_select))

        self._cache.setdefault('brdfile', {})
        if brd_file_path not in self._cache['brdfile']:
            brd_fh = open(brd_file_path, 'rb')
            brd_data = msgpack.unpackb(brd_fh.read())
            brd_fh.close()

            self._cache['brdfile'][brd_file_path] = {}
            self._cache['brdfile'][brd_file_path]['atomgroup'] = bridge.AtomGroup(
                brd_data)

        atomgroup = self._cache['brdfile'][brd_file_path]['atomgroup']
        # selecter = bridge.Select_PathRegex(brd_select)
        selecter = bridge.Select_Path_wildcard(brd_select)
        answer = atomgroup.select(selecter)

        return answer

    def _set_basis_set(self, fragment,
                       basis_set,
                       basis_set_aux=None,
                       basis_set_gridfree=None):
        """
        set basis set to the fragment object.
        """
        if isinstance(fragment, QcFragment):
            for subfrg_name, subfrg in fragment.groups():
                self._set_basis_set(subfrg, basis_set,
                                    basis_set_aux, basis_set_gridfree)
            for atm_name, atm in fragment.atoms():
                atm.basisset = 'O-{}.{}'.format(
                    self._find_basis_set_name(atm, basis_set), atm.symbol)
                if basis_set_aux:
                    atm.basisset_j = 'A-{}.{}'.format(
                        self._find_basis_set_name(atm, basis_set_aux), atm.symbol)
                    atm.basisset_xc = 'A-{}.{}'.format(
                        self._find_basis_set_name(atm, basis_set_aux), atm.symbol)
                else:
                    atm.basisset_j = 'A-{}.{}'.format(
                        self._find_basis_set_name(atm, basis_set), atm.symbol)
                    atm.basisset_xc = 'A-{}.{}'.format(
                        self._find_basis_set_name(atm, basis_set), atm.symbol)
                if basis_set_gridfree:
                    atm.basisset_gridfree = 'O-{}.{}'.format(
                        self._find_basis_set_name(atm, basis_set_gridfree), atm.symbol)
                else:
                    atm.basisset_gridfree = 'O-{}.{}'.format(
                        self._find_basis_set_name(atm, basis_set), atm.symbol)
        else:
            raise QcControlError('Program Error: ', str(fragment))

    def _find_basis_set_name(self, atom, input_obj):
        """
        find suitable basis set name for the atom.

        """
        assert(isinstance(atom, bridge.Atom))
        ans = ''
        if isinstance(input_obj, bridge.basestring):
            ans = bridge.Utils.to_unicode(input_obj)
        elif isinstance(input_obj, dict):
            if atom.name in input_obj:
                ans = input_obj[atom.name]
            elif atom.symbol in input_obj:
                ans = input_obj[atom.symbol]
            else:
                raise QcControlError('type mismatch.', str(input_obj))
        return ans


    # ------------------------------------------------------------------
    # others
    # ------------------------------------------------------------------
    def show_version(self):
        logger.info("=" * 80)
        logger.info("QCLObot version: {version}".format(
            version=str(__version__)))
        logger.info("=" * 80)


if __name__ == '__main__':
    pass
