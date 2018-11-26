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

import os
import sys
import shutil
import logging
from collections import OrderedDict

import proteindf_bridge as bridge
import proteindf_tools as pdf

from .qcatom import QcAtom
from .qcorbitaldata import QcOrbitalData
from .qccommon import get_tmpfile_path
# from .qcframe import QcFrame

logger = logging.getLogger(__name__)

class QcFragment(object):
    # fragmentの密度行列
    _density_matrix_path = 'density.{run_type}.mat'
    # fragmentのLO行列
    _LO_matrix_path = 'LO.{run_type}.mat'
    # fragmentのQCLO行列
    _QCLO_matrix_path = 'QCLO.{run_type}.mat'


    def __init__(self, *args, **kwargs):
        '''
        constructer
        '''
        self._initialize(*args, **kwargs)

        if len(args) > 0:
            rhs = args[0]
            if isinstance(rhs, QcFragment):
                self._copy_constructer(args[0])
            elif isinstance(rhs, bridge.AtomGroup):
                self._construct_by_atomgroup(args[0])
            elif isinstance(rhs, dict):
                self.set_by_raw_data(rhs)

        if 'margin' in kwargs:
            self.margin = kwargs['margin']
        if 'parent' in kwargs:
            self.parent = kwargs['parent']


    def _copy_constructer(self, rhs):
        '''
        copy constructer
        '''
        assert(isinstance(rhs, QcFragment))

        for k, v in rhs.atoms():
            self.set_atom(k, v)
        for k, v in rhs.groups():
            self.set_group(k, v)
        self._name = rhs._name
        self._margin = rhs._margin

        self._parent = rhs._parent
        self._ref_fragment = rhs._ref_fragment


    def _construct_by_atomgroup(self, rhs):
        '''
        '''
        assert(isinstance(rhs, bridge.AtomGroup))

        for k, v in rhs.atoms():
            self.set_atom(k, v)
        for k, v in rhs.groups():
            self.set_group(k, v)
        self._name = rhs.name

        self._ref_fragment = None
        self._parent = None


    def _initialize(self, *args, **kwargs):
        # mandatory variables
        self._atoms = OrderedDict()
        self._groups = OrderedDict()
        self._name = kwargs.get('name', '')
        self._margin = False
        # option variables
        self._ref_fragment = None # isinstance(QcFragment)
        self._parent = None # isinstance(QcFrame, QcFragment)
        # command alias
        self._cmds = self._get_default_cmds()

    def _get_default_cmds(self):
        answer = {}
        answer['mat-extend'] = 'mat-extend'

        return answer

    # state ============================================================
    def get_raw_data(self):
        return self.__get_state__()


    def set_by_raw_data(self, raw_data):
        assert(isinstance(raw_data, dict))
        self.__set_state__(raw_data)


    def __get_state__(self):
        state = {}

        tmp_atoms = []
        for atm_name, atm in self.atoms():
            tmp_atoms.append((atm_name, atm.get_raw_data()))
        state['atoms'] = tmp_atoms

        tmp_grps = []
        for grp_name, grp in self.groups():
            tmp_grps.append((grp_name, grp.get_raw_data()))
        state['groups'] = tmp_grps

        state['name'] = self.name
        state['margin'] = self.margin

        return state


    def __set_state__(self, state):
        assert(isinstance(state, dict))
        self._atoms = OrderedDict()
        if 'atoms' in state:
            for (atm_name, atm_raw) in state.get('atoms'):
                atom = QcAtom()
                atom.set_by_raw_data(atm_raw)
                self.set_atom(atm_name, atom)

        self._groups = OrderedDict()
        if 'groups' in state:
            for (grp_name, grp_raw) in state.get('groups'):
                self.set_group(grp_name, QcFragment(grp_raw, parent=self))

        self._name = state.get('name', '')
        self._margin = state.get('margin', False)


    def _prepare_work_dir(self):
        if not os.path.exists(self.work_dir):
            logger.info('make workdir: {}'.format(self.work_dir))
            os.mkdir(self.work_dir)
        else:
            logger.debug('already exist: {}'.format(self.work_dir))


    # ==================================================================
    # PROPERTIES
    # ==================================================================
    # command alias ----------------------------------------------------
    def set_command_alias(self, cmd_alias_dict):
        for k, v in cmd_alias_dict.items():
            logger.debug("command update: {} -> {}".format(k, v))
            self._cmds[k] = v

    # parent -----------------------------------------------------------
    def _get_parent(self):
        return self._parent
    def _set_parent(self, parent):
        '''
        '''
        from .qcframe import QcFrame
        assert(isinstance(parent, (QcFrame, QcFragment)))
        if (self._parent != None) and (self._parent != parent):
            logger.warn('[{}] parent is overwrite: {} -> {}'.format(
                self.name,
                self._parent.name,
                parent.name))

        self._parent = parent

    parent = property(_get_parent, _set_parent)


    # reference fragment -----------------------------------------------
    def _get_ref_fragment(self):
        return self._ref_fragment
    def _set_ref_fragment(self, frg):
        assert(isinstance(frg, QcFragment))
        logger.info("{header} reference: {ref_parent}/{ref_fragment}".format(
            header=self.header,
            ref_parent=frg.parent.name,
            ref_fragment=frg.name))
        self._ref_fragment = frg
    ref_fragment = property(_get_ref_fragment, _set_ref_fragment)

    # margin ----------------------------------------------------------
    def _get_margin(self):
        return self._margin

    def _set_margin(self, yn):
        self._margin = bool(yn)

    margin = property(_get_margin, _set_margin)


    # work_dir ---------------------------------------------------------
    def _get_work_dir(self):
        '''
        return work_dir path
        '''
        if len(self.name) == 0:
            logger.critical('fragment.name is not define: {}'.format(repr(self.name)))
            raise

        parent_path = ''
        if self.parent != None:
            parent_path = self.parent.work_dir
        else:
            logger.debug('not set parent.')

        return os.path.join(parent_path, self.name)

    work_dir = property(_get_work_dir)

    def _check_path(self, path):
        if not os.path.exists(path):
            logger.warn("{header} NOT FOUND: {path}".format(
                header=self.header, path=path))

    def _check_matrix(self, path, expect_row, expect_col):
        answer = False
        if pdf.Matrix.is_loadable(path):
            (row, col) = pdf.Matrix.get_size(path)
            logger.warning("actual({}, {}) <> expect({}, {})".format(row, col, expect_row, expect_col))
            if (((expect_row == None) or (expect_row == row)) and
                ((expect_col == None) or (expect_col == col))):
                answer = True
        return answer

    def _check_symmetric_matrix(self, path, dim):
        answer = False
        if pdf.SymmetricMatrix.is_loadable(path):
            (row, col) = pdf.SymmetricMatrix.get_size(path)
            if row == dim and col == dim:
                answer = True
        return answer

    # ==================================================================
    # molecular properties
    # ==================================================================
    # name -------------------------------------------------------------
    def _get_name(self):
        return self._name
    def _set_name(self, name):
        name = bridge.Utils.to_unicode(name)
        self._name = name
    name = property(_get_name, _set_name)


    # number of AOs ----------------------------------------------------
    def get_number_of_AOs(self):
        def get_number_of_AOs_sub(ag):
            AOs = 0
            for key, subgrp in ag.groups():
                AOs += get_number_of_AOs_sub(subgrp)
            for key, atm in ag.atoms():
                AOs += atm.get_number_of_AOs()
            return AOs
        return get_number_of_AOs_sub(self)


    # AtomGroup -------------------------------------------------------
    def get_AtomGroup(self):
        ag = bridge.AtomGroup()
        for subgrp_name, subgrp in self.groups():
            ag.set_group(subgrp_name, subgrp.get_AtomGroup())
        for atm_name, atm in self.atoms():
            a = atm.get_Atom()
            ag.set_atom(atm_name, atm)
        return ag


    # basisset ---------------------------------------------------------
    def set_basisset(self, pdfparam):
        for key, frg in self.groups():
            frg.set_basisset(pdfparam)
        for key, atm in self.atoms():
            symbol = atm.symbol
            if symbol != 'X':
                atomlabel = atm.atomlabel
                bsname = atm.basisset
                bsname_j = atm.basisset_j
                bsname_xc = atm.basisset_xc
                bsname_gridfree = atm.basisset_gridfree
                pdfparam.set_basisset_name(atomlabel, bsname)
                pdfparam.set_basisset_j_name(atomlabel, bsname_j)
                pdfparam.set_basisset_xc_name(atomlabel, bsname_xc)
                pdfparam.set_basisset_gridfree_name(atomlabel, bsname_gridfree)


    # orbital info -----------------------------------------------------
    def get_orbital_info(self):
        basis2 = pdf.Basis2()

        orbital_info = []
        for subgrp_key, subgrp in self.groups():
            orbital_info.extend(subgrp.get_orbital_info())

        for atom_key, qc_atom in self.atoms():
            atom = QcAtom(qc_atom)
            basisset_name = qc_atom.basisset
            basisset = basis2.get_basisset(basisset_name)
            num_of_CGTOs = len(basisset)
            for CGTO_index in range(num_of_CGTOs):
                CGTO = basisset[CGTO_index]
                shell_type = CGTO.shell_type
                shell_type_id = pdf.ContractedGTO.get_shell_type_id(shell_type)
                num_of_basis_type = shell_type_id * 2 + 1
                for basis_type in range(num_of_basis_type):
                    data = QcOrbitalData(atom = atom,
                                         basisset_name = basisset_name,
                                         CGTO_index = CGTO_index,
                                         basis_type = basis_type)
                    orbital_info.append(data)
        return orbital_info


    # ==================================================================
    # modeling
    # ==================================================================
    # ==================================================================
    # atom
    # ==================================================================
    def get_number_of_atoms(self):
        return len(self._atoms)


    def get_number_of_all_atoms(self):
        answer = 0
        for key, grp in self.groups():
            answer += grp.get_number_of_all_atoms()
        answer += len(self._atoms)
        return answer


    def sum_of_atomic_number(self):
        """
        原子数の総和を返す
        """
        answer = 0.0
        for key, grp in self.groups():
            answer += grp.sum_of_atomic_number()
        for key, atm in self.atoms():
            answer += atm.atomic_number
        return answer


    def get_atom_kinds(self):
        """
        原子種(シンボル)のリストを返す
        """
        answer = set()
        for key, group in self.groups():
            tmp = group.get_atom_kinds()
            answer.update(tmp)
        for key, atom in self.atoms():
            answer.add(atom.symbol)
        return answer


    def get_atom(self, key_or_name):
        '''
        入力されたkeyもしくは名前の原子が含まれている場合、その原子を返す。
        無い場合はNoneを返す。
        '''
        if key_or_name in self._atoms:
            return self._atoms.get(key_or_name, None)
        else:
            for k, atm in self.atoms():
                if atm.name == key_or_name:
                    return atm
        return None


    def set_atom(self, key, value):
        '''
        Set QcAtom object.

        '''
        self._atoms[key] = QcAtom(value)


    def atoms(self):
        '''
        原子のリストを返す
        '''
        for k, v in self._atoms.items(): # based on collections.OrderedDict
            yield(k, v)

    # ==================================================================
    # group
    # ==================================================================
    def get_number_of_groups(self):
        return len(self._groups)


    def set_group(self, key, fragment):
        '''
        Set QcFragment object.

        This method override the bridge.AtomGroup method.
        '''
        self._groups[key] = QcFragment(fragment)
        self._groups[key].parent = self


    def groups(self):
        '''
        '''
        for k, v in self._groups.items(): # based on collections.OrderedDict
            yield(k, v)


    # def _delete_group(self, key):
    #    self._groups.pop(key)


    def grouping_subfragments(self):
        '''
        子グループを自分の原子リストに組入れ、その子グループを削除する
        '''
        # set atoms in subgroup to my atom list
        for key_subgrp, subgrp in self.groups():
            subgrp.grouping_subfragments()
            for key_atom, atom in subgrp.atoms():
                new_key_atom = '{}/{}'.format(key_subgrp, key_atom)
                self.set_atom(new_key_atom, atom)
        # delete subgroup
        self._groups = OrderedDict()


    # ==================================================================
    # guess density
    # ==================================================================
    # set corresponding density matrix ---------------------------------
    def set_density_matrix(self, in_path, run_type='rks'):
        abs_in_path = os.path.abspath(in_path)
        density_matrix_path = os.path.abspath(self._get_density_matrix_path(run_type))
        if abs_in_path != density_matrix_path:
            self._prepare_work_dir()
            logger.info("{header} save density matrix as {path}".format(header=self.header,
                                                                        path=density_matrix_path))
            shutil.move(abs_in_path,
                        density_matrix_path)
        else:
            logger.warning('not set the same density matrix')


    def get_density_matrix(self, run_type="rks"):
        """本フラグメントの密度行列ファイルのパスを返す。

        ファイルが存在しない場合は親フレーム分子に作成を依頼する。
        """
        path = self._get_density_matrix_path(run_type)
        if os.path.isfile(path) != True:
            self.parent.pickup_density_matrix(run_type)
        assert(os.path.isfile(path))
        return path


    def _get_density_matrix_path(self, run_type):
        """本フラグメントの密度行列ファイルのパスを返す
        """
        return os.path.join(self.work_dir, self._density_matrix_path.format(run_type=run_type))


    # create density matrix ---------------------------------------------
    def prepare_guess_density_matrix(self, run_type):
        '''
        guess_densityに必要な密度行列のパスを返す
        subgroupを持っている場合はマージした密度行列を作成し、そのパスを返す
        '''
        self._prepare_work_dir()

        logger.info("{header} prepare guess density matrix: start".format(header=self.header))
        guess_density_matrix_path = os.path.join(self.work_dir,
                                                 'guess.density.{}.{}.mat'.format(run_type, self.name))

        # 既存のデータを消去する
        if os.path.exists(guess_density_matrix_path):
            logger.info("{header}/remove {path}".format(
                header=self.header, path=guess_density_matrix_path))
            os.remove(guess_density_matrix_path)
        # create new matrix file
        mat = pdf.SymmetricMatrix()
        mat.save(guess_density_matrix_path)

        # subgroup
        logger.info('{header} get subgrp density matrix'.format(header=self.header))
        num_of_AOs_subgrp = 0
        for subgrp_name, subgrp in self.groups():
            logger.info('{header} subgroup name={subgrp_name}'.format(
                header=self.header, subgrp_name=subgrp_name))

            subgrp.set_command_alias(self._cmds)
            subgrp_guess_density_matrix_path = subgrp.prepare_guess_density_matrix(run_type)
            if not os.path.exists(subgrp_guess_density_matrix_path):
                logger.warn('NOT found: subgrp.guess.dens.mat={}'.format(subgrp_guess_density_matrix_path))
                continue
            self._check_path(subgrp_guess_density_matrix_path)
            assert(self._check_symmetric_matrix(guess_density_matrix_path, num_of_AOs_subgrp))

            # merge subgrp to main density matrix
            subgrp_AOs = subgrp.get_number_of_AOs()
            assert(self._check_symmetric_matrix(subgrp_guess_density_matrix_path, subgrp_AOs))
            num_of_AOs_subgrp += subgrp_AOs

            logger.debug('(sub) {} -d '.format(self._cmds['mat-extend']))
            logger.debug('    {}'.format(guess_density_matrix_path))
            logger.debug('    {}'.format(subgrp_guess_density_matrix_path))
            logger.debug('    {}'.format(guess_density_matrix_path))
            pdf.run_pdf([self._cmds['mat-extend'], '-d',
                         guess_density_matrix_path,
                         subgrp_guess_density_matrix_path,
                         guess_density_matrix_path])
        assert(self._check_symmetric_matrix(guess_density_matrix_path, num_of_AOs_subgrp))

        # self
        logger.info('{header} get self density matrix'.format(header=self.header))
        if len(self._atoms) > 0:
            # (計算済みの)参照元の密度行列パスを取得する
            #  parentは未計算(これから計算)なので密度行列は取得できない。
            assert(self.ref_fragment != None)
            my_density_matrix_path = self.ref_fragment.get_density_matrix(run_type)
            logger.info('{header} reference density matrix path: {path}'.format(
                header=self.header, path=my_density_matrix_path))
            self._check_path(my_density_matrix_path)

            logger.debug('{} -d {} {} {}'.format(
                self._cmds['mat-extend'],
                guess_density_matrix_path, my_density_matrix_path, guess_density_matrix_path))
            pdf.run_pdf([self._cmds['mat-extend'], '-d',
                         guess_density_matrix_path,
                         my_density_matrix_path,
                         guess_density_matrix_path])

        # check
        self._check_path(guess_density_matrix_path)
        assert(self._check_symmetric_matrix(guess_density_matrix_path, self.get_number_of_AOs()))

        logger.info("{header} prepare guess density matrix: end".format(header=self.header))
        return guess_density_matrix_path



    # LO --------------------------------------------------------------
    def set_LO_matrix(self, in_path, run_type='rks'):
        abs_in_path = os.path.abspath(in_path)
        LO_matrix_path = os.path.abspath(self.get_LO_matrix_path(run_type))
        if abs_in_path != LO_matrix_path:
            self._prepare_work_dir()
            logger.info("{header} save LO matrix as {path}".format(header=self.header,
                                                                   path=LO_matrix_path))
            shutil.move(abs_in_path,
                        LO_matrix_path)
        else:
            logger.warning('not set the same LO matrix')

    def get_LO_matrix_path(self, run_type):
        return os.path.join(self.work_dir, self._LO_matrix_path.format(run_type=run_type))


    # QCLO ------------------------------------------------------------
    def set_QCLO_matrix(self, in_path, run_type='rks'):
        abs_in_path = os.path.abspath(in_path)
        QCLO_matrix_path = os.path.abspath(self._get_QCLO_matrix_path(run_type))
        if abs_in_path != QCLO_matrix_path:
            self._prepare_work_dir()
            logger.info("{header} save QCLO matrix as {path}".format(header=self.header,
                                                                     path=QCLO_matrix_path))
            shutil.move(abs_in_path,
                        QCLO_matrix_path)
        else:
            logger.warning('not set the same QCLO matrix')


    def get_QCLO_matrix_path(self, run_type="rks", force=False):
        """本フラグメントのQCLO行列ファイルのパスを返す。

        ファイルが存在しない場合は親フレーム分子に作成を依頼する。
        """
        path = self._get_QCLO_matrix_path(run_type)
        if os.path.isfile(path) != True:
            self.parent.pickup_QCLO_matrix(run_type, force)
        assert(os.path.isfile(path))
        return path


    def _get_QCLO_matrix_path(self, run_type):
        """本フラグメントのQCLO行列ファイルのパスを返す
        """
        return os.path.join(self.work_dir, self._QCLO_matrix_path.format(run_type=run_type))


    def prepare_guess_QCLO_matrix(self, run_type, request_frame, force=False):
        '''
        prepare QCLO matrix

        QCLO行列のパスを返す
        subgroupを持っている場合はマージしたQCLO行列を作成し、そのパスを返す
        '''
        self._prepare_work_dir()

        logger.info("{header} prepare QCLO matrix: start".format(header=self.header))
        guess_QCLO_matrix_path = get_tmpfile_path()

        request_orbinfo = request_frame.get_orbital_info()
        request_num_of_AOs = request_frame.get_number_of_AOs()

        # 既存のデータを消去する
        if os.path.isfile(guess_QCLO_matrix_path):
            logger.warning("{header} remove existed fragment QCLO file: {path}".format(
                header=self.header, path=guess_QCLO_matrix_path))
            os.remove(guess_QCLO_matrix_path)

        # subgroup
        logger.info("{header} get subgroup QCLO matrix".format(header=self.header))
        for subgrp_name, subgrp in self.groups():
            logger.info('{header} subgroup name={subgrp_name}'.format(
                header=self.header, subgrp_name=subgrp_name))

            subgrp.set_command_alias(self._cmds)
            subgrp_guess_QCLO_matrix_path = subgrp.prepare_guess_QCLO_matrix(run_type, request_frame)
            if not os.path.exists(subgrp_guess_QCLO_matrix_path):
                logger.warn('NOT found: subgrp.guess.QCLO.mat={}'.format(subgrp_guess_QCLO_matrix_path))
                continue
            self._check_path(subgrp_guess_QCLO_matrix_path)

            # 行数は変えずに列方向に追加("pdf-mat-ext -c")
            logger.debug('{} -c '.format(self._cmds['mat-extend']))
            logger.debug('    {}'.format(guess_QCLO_matrix_path))
            logger.debug('    {}'.format(subgrp_guess_QCLO_matrix_path))
            logger.debug('    {}'.format(guess_QCLO_matrix_path))
            pdf.run_pdf([self._cmds['mat-extend'], '-c',
                         guess_QCLO_matrix_path,
                         subgrp_guess_QCLO_matrix_path,
                         guess_QCLO_matrix_path])
            self._check_path(guess_QCLO_matrix_path)
            assert(self._check_matrix(guess_QCLO_matrix_path, request_num_of_AOs, None))

        # 自分のQCLO情報
        if len(self._atoms) > 0:
            logger.info("{header} get self QCLO matrix".format(header=self.header))
            ref_frame = self.ref_fragment.parent
            ref_orbinfo = ref_frame.get_orbital_info()
            ref_num_of_AOs = ref_frame.get_number_of_AOs()

            my_qclo_matrix_path = self.ref_fragment.get_QCLO_matrix_path(run_type, force)
            logger.info('{header} reference QCLO matrix path: {path}'.format(
                header=self.header, path=my_qclo_matrix_path))
            self._check_path(my_qclo_matrix_path)

            QCLO_mat = pdf.Matrix()
            QCLO_mat.load(my_qclo_matrix_path)
            if QCLO_mat.rows != parent_num_of_AOs:
                logger.warning("QCLO matrix row(= {qclo_row}) is not equal to the parent AOs(= {ao})".format(
                    qclo_row = QCLO_mat.row,
                    ao = parent_num_of_AOs))
            num_of_MOs = QCLO_mat.cols
            guess_QCLO_mat = pdf.Matrix(request_num_of_AOs, num_of_MOs)

            for request_AO_index in range(request_num_of_AOs):
                for ref_AO_index in range(ref_num_of_AOs):
                    if request_orbinfo[request_AO_index] == ref_orbinfo[ref_AO_index]:
                        for MO_index in range(num_of_MOs):
                            v = QCLO_mat.get(ref_AO_index, MO_index)
                            guess_QCLO_mat.set(request_AO_index, MO_index, v)
            my_guess_QCLO_matrix_path = os.path.join(self.work_dir, "guess_QCLO.part.mat")
            guess_QCLO_mat.save(my_guess_QCLO_matrix_path)

            # 行数は変えずに列方向に追加("pdf-mat-ext -c")
            logger.debug('{} -c '.format('mat-extend'))
            logger.debug('    {}'.format(guess_QCLO_matrix_path))
            logger.debug('    {}'.format(my_guess_QCLO_matrix_path))
            logger.debug('    {}'.format(guess_QCLO_matrix_path))
            pdf.run_pdf([self._cmds['mat-extend'], '-c',
                         guess_QCLO_matrix_path,
                         my_guess_QCLO_matrix_path,
                         guess_QCLO_matrix_path])
            self._check_path(guess_QCLO_matrix_path)
        else:
            logger.info("{header} no belonging atoms found. No QCLO created.")

        # check
        self._check_path(guess_QCLO_matrix_path)
        assert(self._check_matrix(guess_QCLO_matrix_path, request_num_of_AOs, None))

        logger.info("{header} prepare QCLO: end".format(header=self.header))
        return guess_QCLO_matrix_path


    # ==================================================================
    # operator
    # ==================================================================
    def __getitem__(self, key):
        '''
        operator[] for getter

        keyが一致した原子団、原子を返す。
        もしkeyが一致しなければ、名前から検索する。
        '''
        if (self.has_group(key) == True):
            return self._groups[key]
        elif key in self._atoms:
            return self._atoms[key]
        else:
            for k, grp in self.groups():
                if grp.name == key:
                    return grp
            for k, atm in self.atoms():
                if atm.name == key:
                    return atm
        raise KeyError(key)


    def __setitem__(self, key, value):
        '''
        operator[] for setter
        '''
        if (isinstance(value, (QcFragment, bridge.AtomGroup)) == True):
            self.set_group(key, value)
        elif (isinstance(value, (QcAtom, bridge.Atom)) == True):
            self.set_atom(key, value)
        else:
            raise ValueError(value)


    # operator == ------------------------------------------------------
    def __eq__(self, rhs):
        if (rhs == None) or (isinstance(rhs, QcFragment) == False):
            return False
        return ((self.parent == rhs.parent) and
                (self.name == rhs.name))


    def __ne__(self, rhs):
        return not self.__eq__(rhs)


    # str -------------------------------------------------------------
    def __str__(self):
        answer = ''
        for key, subgrp in self.groups():
            answer += '>>>> {}\n'.format(key)
            answer += str(subgrp) + '\n'
        for key, atom in self.atoms():
            answer += 'k:{} {}\n'.format(key, str(atom))

        return answer


    # ==================================================================
    # debug
    # ==================================================================
    def _get_logger_header(self):
        header = ""
        if self.parent != None:
            header += "{}".format(self.parent.name)
        header += "/{}>".format(self.name)
        return header
    header=property(_get_logger_header)
