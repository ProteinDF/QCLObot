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

import pdfbridge as bridge
import pdfpytools as pdf
import qclobot as qclo

class QcFragment(object):
    def __init__(self, *args, **kwargs):
        '''
        constructer
        '''
        self._logger = logging.getLogger(__name__)

        self._atoms = OrderedDict()
        self._groups = OrderedDict()
        self._name = kwargs.get('name', '')
        self._qc_parent = None
        self._density_matrix_path = None
        self._Clo_path = None
        self._QCLO_matrix_path = None

        if len(args) > 0:
            rhs = args[0]
            if isinstance(rhs, QcFragment):
                self._copy_constructer(args[0])
            elif isinstance(rhs, bridge.AtomGroup):
                self._transform_atomgroup(args[0])

        if 'qc_parent' in kwargs:
            self.qc_parent = kwargs['qc_parent']
                    
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
        self._qc_parent = rhs._qc_parent
        self._density_matrix_path = rhs._density_matrix_path
        self._Clo_path = rhs._Clo_path
        self._QCLO_matrix_path = rhs._QCLO_matrix_path

    def _transform_atomgroup(self, rhs):
        '''
        '''
        assert(isinstance(rhs, bridge.AtomGroup))

        for k, v in rhs.atoms():
            self.set_atom(k, v)
        for k, v in rhs.groups():
            self.set_group(k, v)
        self._name = rhs.name
        self._qc_parent = None
        self._density_matrix_path = None
        self._Clo_path = None
        self._QCLO_matrix_path = None
            
        
    def _prepare_work_dir(self):
        if not os.path.exists(self.work_dir):
            self._logger.info('make workdir: {}'.format(self.work_dir))
            os.mkdir(self.work_dir)
        else:
            self._logger.debug('already exist: {}'.format(self.work_dir))
        
    # ==================================================================
    # PROPERTIES
    # ==================================================================
    # QcParent ---------------------------------------------------------
    def _get_qc_parent(self):
        return self._qc_parent

    def _set_qc_parent(self, qc_parent):
        '''
        '''
        from qcframe import QcFrame
        assert(isinstance(qc_parent, (QcFrame, QcFragment)))
        if (self.qc_parent != None) and (self.qc_parent != qc_parent):
            self._logger.warn('[{}] qc_parent is overwrite: {} -> {}'.format(
                self.name,
                self._qc_parent.name,
                qc_parent.name))
        self._qc_parent = qc_parent

    qc_parent = property(_get_qc_parent, _set_qc_parent)
        
    # work_dir ---------------------------------------------------------
    def _get_work_dir(self):
        '''
        return work_dir path
        '''
        assert(len(self.name) > 0)

        parent_path = ''
        if self.qc_parent != None:
            parent_path = self.qc_parent.work_dir
        else:
            self._logger.debug('not set parent.')

        return os.path.join(parent_path, self.name)

    work_dir = property(_get_work_dir)

    def _check_path(self, path):
        if not os.path.exists(path):
            self._logger.warn('NOT FOUND: {}'.format(path))
        
    # ==================================================================
    # molecular properties
    # ==================================================================
    # name -------------------------------------------------------------
    def _get_name(self):
        return self._name

    def _set_name(self, name):
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
            atom = qclo.QcAtom(qc_atom)
            basisset_name = qc_atom.basisset
            basisset = basis2.get_basisset(basisset_name)
            num_of_CGTOs = len(basisset)
            for CGTO_index in range(num_of_CGTOs):
                CGTO = basisset[CGTO_index]
                shell_type = CGTO.shell_type
                shell_type_id = pdf.ContractedGTO.get_shell_type_id(shell_type)
                num_of_basis_type = shell_type_id * 2 + 1
                for basis_type in range(num_of_basis_type):
                    data = qclo.QcOrbitalData(atom = atom,
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
        print('set_atom: {}: {}'.format(key, str(value)))
        self._atoms[key] = qclo.QcAtom(value, qc_parent = self)

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
        if fragment.qc_parent == None:
            self._groups[key].qc_parent = self

    def groups(self):
        '''
        '''
        for k, v in self._groups.items(): # based on collections.OrderedDict
            yield(k, v)

    # ==================================================================
    # guess density
    # ==================================================================
    # set corresponding density matrix ---------------------------------
    def set_density_matrix(self, in_path):
        basename = os.path.basename(in_path)

        input_abs_path = os.path.abspath(in_path)
        my_abs_path = os.path.abspath(os.path.join(self.work_dir, basename))
        if input_abs_path != my_abs_path:
            self._prepare_work_dir()
            shutil.copy(os.path.abspath(in_path),
                        os.path.abspath(my_abs_path))
        self._density_matrix_path = basename

        self._logger.info('set dens. mat. path: [{name}]/{path}'.format(
            name = self.name,
            path = os.path.join(self.work_dir,
                                self._density_matrix_path)))
        
    # create density matrix ---------------------------------------------
    def get_guess_density_matrix(self, run_type):
        '''
        guess_densityに必要な密度行列のパスを返す
        subgroupを持っている場合はマージした密度行列を作成し、そのパスを返す
        '''
        self._logger.info('>>>> get_guess_density_matrix: {}/{}'.format(self.qc_parent.name, self.name))
        guess_density_matrix_path = os.path.join(self.work_dir,
                                                 'guess.density.{}.{}.mat'.format(run_type, self.name))
        
        # サブユニットからの行列をこれから追記するので
        # 既存のデータを消去する
        if os.path.exists(guess_density_matrix_path):
            os.remove(guess_density_matrix_path)
        
        for subgrp_name, subgrp in self.groups():
            subgrp_guess_density_matrix_path = subgrp.get_guess_density_matrix(run_type)
            if not os.path.exists(subgrp_guess_density_matrix_path):
                self._logger.warn('NOT found: subgrp.guess.dens.mat={}'.format(subgrp_guess_density_matrix_path))
                continue
            self._check_path(subgrp_guess_density_matrix_path)
            
            # merge subgrp to main density matrix
            (is_loadable, row1, col1) = pdf.SymmetricMatrix.is_loadable(guess_density_matrix_path)
            assert(is_loadable == True)
            (is_loadable, row2, col2) = pdf.SymmetricMatrix.is_loadable(subgrp_guess_density_matrix_path)
            assert(is_loadable == True)

            self._logger.info('(sub) mat-ext -d ')
            self._logger.info('    {}'.format(guess_density_matrix_path))
            self._logger.info('    {}'.format(subgrp_guess_matrix_path))
            self._logger.info('    {}'.format(guess_density_matrix_path))
            pdf.run_pdf(['mat-ext', '-d',
                         guess_density_matrix_path,
                         subgrp_guess_density_matrix_path,
                         guess_density_matrix_path])
            self._check_path(guess_density_matrix_path)
            (is_loadable, row3, col3) = pdf.SymmetricMatrix.is_loadable(guess_density_matrix_path)
            assert(is_loadable == True)
            assert(row1 + row2 == row3)
            assert(col1 + col2 == col3)

        if self._density_matrix_path != None:
            path = os.path.join(self.work_dir, self._density_matrix_path)
            self._check_path(path)
            self._logger.info('mat-ext -d ')
            self._logger.info('    {}'.format(guess_density_matrix_path))
            self._logger.info('    {}'.format(path))
            self._logger.info('    {}'.format(guess_density_matrix_path))
            pdf.run_pdf(['mat-ext', '-d',
                         guess_density_matrix_path,
                         path,
                         guess_density_matrix_path])

        self._check_path(guess_density_matrix_path)
        self._logger.info('<<<< get_guess_density_matrix: {}/{}'.format(self.qc_parent.name, self.name))
        return guess_density_matrix_path


    def get_guess_density_matrix0(self, run_type):
        '''
        guess_densityに必要な密度行列のパスを返す
        subgroupを持っている場合はマージした密度行列を作成し、そのパスを返す
        '''
        guess_matrix_path = os.path.join(self.work_dir,
                                         'guess.Ppq.{}.{}.mat'.format(run_type, self.name))
        
        # サブユニットからの行列をこれから追記するので
        # 既存のデータを消去する
        if os.path.exists(guess_matrix_path):
            os.remove(guess_matrix_path)
        
        for subgrp_name, subgrp in self.groups():
            self._logger.info('{} merge {} from the subgroup, {}'.format(
                __name__,
                guess_matrix_path,
                subgrp_guess_matrix_path))

            subgrp_guess_matrix_path = subgrp.get_guess_density_matrix(run_type)

            (is_loadable, row1, col1) = pdf.SymmetricMatrix.is_loadable(density_matrix_path)
            assert(is_loadable == True)
            (is_loadable, row2, col2) = pdf.SymmetricMatrix.is_loadable(subgrp_density_matrix_path)
            assert(is_loadable == True)

            self._logger.info('mat-ext -d ')
            self._logger.info('    {}'.format(guess_matrix_path))
            self._logger.info('    {}'.format(subgrp_guess_matrix_path))
            self._logger.info('    {}'.format(guess_matrix_path))
            pdf.run_pdf(['mat-ext', '-d',
                         guess_matrix_path,
                         subgrp_guess_matrix_path,
                         guess_matrix_path])
            (is_loadable, row3, col3) = pdf.SymmetricMatrix.is_loadable(density_matrix_path)
            assert(is_loadable == True)
            assert(row1 + row2 == row3)
            assert(col1 + col2 == col3)

        if self._density_matrix_path != None:
            path = os.path.join(self.work_dir, self._density_matrix_path)
            #self._logger.info('{} ext(self): {} copy from {}'.format(
            #    __name__,
            #    guess_matrix_path,
            #    path))
            self._logger.info('mat-ext -d ')
            self._logger.info('    {}'.format(guess_matrix_path))
            self._logger.info('    {}'.format(path))
            self._logger.info('    {}'.format(guess_matrix_path))
            pdf.run_pdf(['mat-ext', '-d',
                         guess_matrix_path,
                         path,
                         guess_matrix_path])

        return guess_matrix_path

    # LO --------------------------------------------------------------
    def _get_Clo_path(self):
        path = self._Clo_path
        if path != None:
            if ((self.qc_parent) and (self.qc_parent.work_dir)):
                path = os.path.join(self.qc_parent.work_dir, path)
        return path

    def _set_Clo_path(self, path):
        self._Clo_path = str(path)

    Clo_path = property(_get_Clo_path, _set_Clo_path)
        
    # QCLO ------------------------------------------------------------
    def set_QCLO_matrix(self, in_path):
        basename = os.path.basename(in_path)

        input_abs_path = os.path.abspath(in_path)
        my_abs_path = os.path.abspath(os.path.join(self.work_dir, basename))
        if input_abs_path != my_abs_path:
            self._prepare_work_dir()
            shutil.copy(os.path.abspath(in_path),
                        os.path.abspath(my_abs_path))
        self._QCLO_matrix_path = basename

        self._logger.info('set QCLO mat. path: [{name}]/{path}'.format(
            name = self.name,
            path = os.path.join(self.work_dir,
                                self._QCLO_matrix_path)))
        
    def _get_QCLO_matrix_path(self):
        self._logger.info('_get_QCLO_matrix_path() ')
        self._logger.info('  work_dir: {}'.format(self.work_dir))
        self._logger.info('  _QCLO_matrix_path: {}'.format(self._QCLO_matrix_path))

        answer = None
        if self._QCLO_matrix_path != None:
            answer = os.path.abspath(os.path.join(self.work_dir, self._QCLO_matrix_path))
        return answer
    
    #def _set_QCLO_path(self, path):
    #    self._QCLO_path = str(path)
    #
    QCLO_matrix_path = property(_get_QCLO_matrix_path)

    def prepare_QCLO_matrix(self, run_type, request_frame):
        '''
        prepare QCLO matrix
        '''
        self._prepare_work_dir()
        guess_QCLO_matrix_path = os.path.join(self.work_dir, 'guess_QCLO.mat')

        self._logger.info('>>>> prepare QCLO: {}/{}'.format(self.qc_parent.name, self.name))
        self._logger.info('     for {}'.format(guess_QCLO_matrix_path))

        # サブユニットからの行列をこれから追記するので
        # 既存のデータを消去する
        if os.path.exists(guess_QCLO_matrix_path):
            os.remove(guess_QCLO_matrix_path)

        request_orbinfo = request_frame.get_orbital_info()

        # subgroup
        for subgrp_name, subgrp in self.groups():
            subgrp_guess_QCLO_matrix_path = subgrp.prepare_QCLO_matrix(run_type, request_frame)
            
            self._check_path(subgrp_guess_QCLO_matrix_path)
            self._logger.info('mat-ext -c ')
            self._logger.info('    {}'.format(guess_QCLO_matrix_path))
            self._logger.info('    {}'.format(subgrp_guess_QCLO_matrix_path))
            self._logger.info('    {}'.format(guess_QCLO_matrix_path))
            pdf.run_pdf(['mat-ext', '-c',
                         guess_QCLO_matrix_path,
                         subgrp_guess_QCLO_matrix_path,
                         guess_QCLO_matrix_path])
            self._check_path(guess_QCLO_matrix_path)
            self._logger.info('')
            
            
        # 自分のQCLO情報
        parent_frame = self.qc_parent
        parent_orbinfo = parent_frame.get_orbital_info()
        request_num_of_AOs = request_frame.get_number_of_AOs()
        parent_num_of_AOs = parent_frame.get_number_of_AOs()

        self._logger.info('QCLO_matrix_path: {}, parent={}'.format(self.QCLO_matrix_path,
                                                                   self.qc_parent.name))
        if self.QCLO_matrix_path != None:
            QCLO_mat = pdf.Matrix()

            self._check_path(self.QCLO_matrix_path)

            QCLO_mat.load(self.QCLO_matrix_path)
            assert(QCLO_mat.rows == parent_num_of_AOs)
            num_of_MOs = QCLO_mat.cols
            guess_QCLO_mat = pdf.Matrix(request_num_of_AOs,
                                        num_of_MOs)
        
            for request_AO_index in range(request_num_of_AOs):
                for parent_AO_index in range(parent_num_of_AOs):
                    if request_orbinfo[request_AO_index] == parent_orbinfo[parent_AO_index]:
                        for MO_index in range(num_of_MOs):
                            v = QCLO_mat.get(parent_AO_index, MO_index)
                            guess_QCLO_mat.set(request_AO_index, MO_index, v)
            my_guess_QCLO_matrix_path = os.path.join(self.work_dir,
                                                     'guess_QCLO.part.mat')
            guess_QCLO_mat.save(my_guess_QCLO_matrix_path)

            self._logger.info('mat-ext -c ')
            self._logger.info('    {}'.format(guess_QCLO_matrix_path))
            self._logger.info('    {}'.format(my_guess_QCLO_matrix_path))
            self._logger.info('    {}'.format(guess_QCLO_matrix_path))
            pdf.run_pdf(['mat-ext', '-c',
                         guess_QCLO_matrix_path,
                         my_guess_QCLO_matrix_path,
                         guess_QCLO_matrix_path])
            self._check_path(guess_QCLO_matrix_path)
            self._logger.info('')

        self._logger.info('prepare QCLO: name={}: {}'.format(self.name,
                                                             guess_QCLO_matrix_path))

        self._check_path(guess_QCLO_matrix_path)
        self._logger.info('<<<< prepare QCLO: {}/{}\n'.format(self.qc_parent.name, self.name))
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
        elif (isinstance(value, (qclo.QcAtom, bridge.Atom)) == True):
            self.set_atom(key, value)
        else:
            raise ValueError(value)


    # operator == ------------------------------------------------------
    def __eq__(self, rhs):
        if rhs == None:
            return False
        return ((self.qc_parent == rhs.qc_parent) and
                (self.name == rhs.name))
        
    def __ne__(self, rhs):
        return not self.__eq__(rhs)

    # str -------------------------------------------------------------
    def __str__(self):
        answer = ''
        for key, subgrp in self.groups():
            answer += '>>>> {}'.format(key)
            answer += str(subgrp) + '\n'
        for key, atom in self.atoms():
            answer += 'k:{} {}\n'.format(key, str(atom))

        return answer 
    

