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
import logging
import math
from collections import OrderedDict
import shutil

try:
    import msgpack
except:
    import msgpack_pure as msgpack

import pdfbridge as bridge
import pdfpytools as pdf
import qclobot as qclo

logger = logging.getLogger(__name__)

class QcFrame(object):
    _pdfparam_filename = 'pdfparam.mpac'
    _db_filename = 'pdfresults.db'

    # ------------------------------------------------------------------
    def __init__(self, name, *args, **kwargs):
        '''
        create QcFrame object.
        
        name: name of the frame molecule
        '''
        if kwargs.get('debug'):
            logger.addHandler(logging.StreamHandler())
            logger.setLevel(logging.DEBUG)
        else:
            logger.addHandler(logging.NullHandler())
            logger.setLevel(logging.INFO)

        # mandatory parameter
        self._name = name
        self._fragments = OrderedDict()
        self._charge = 0
        self._state = {} # 状態保持
        
        self._prepare_work_dir()
        self._load()

        # cache data
        self._cache = {}

        if ((len(args) > 0) and isinstance(args[0], QcFrame)):
            self._copy_constructer(args[0])
            
    # copy constructer
    def _copy_constructer(self, rhs):
        self._name = rhs._name
        self._fragments = copy.deepcopy(rhs._fragments)
        self._charge = rhs._charge
        self._state = copy.deepcopy(rhs._state)

    def __del__(self):
        self._save()

    # save & load ------------------------------------------------------
    def _load(self):
        path = os.path.join(self.work_dir, 'qcframe.mpac')
        if os.path.exists(path):
            logger.info('load the fragment state: {}'.format(path))
            f = open(path, 'rb')
            packed = f.read()
            state_dat = msgpack.unpackb(packed)
            f.close()
            state_dat = bridge.Utils.to_unicode_dict(state_dat)
            self.set_by_raw_data(state_dat)
        else:
            logger.debug('not found the state file')
        
    def _save(self):
        path = os.path.join(self.work_dir, 'qcframe.mpac')
        # logger.info('save the fragment state: {}'.format(path))

        state_dat = self.get_raw_data()
        packed = msgpack.packb(state_dat)
        f = open(path, 'wb')
        f.write(packed)
        f.close()

    def get_raw_data(self):
        return self.__getstate__()
        
    def set_by_raw_data(self, raw_data):
        self.__setstate__(raw_data)
        
    def __getstate__(self):
        state = {}
        state['name'] = self.name
        tmp_frgs = []
        for k, frg in self.fragments():
            tmp_frgs.append((k, frg.get_raw_data()))
        state['fragments'] = tmp_frgs
        state['charge'] = self.charge
        state['state'] = self._state
        return state

    def __setstate__(self, state):
        assert(isinstance(state, dict))
        self._name = state.get('name')
        self._fragments = OrderedDict()
        if 'fragments' in state:
            for (k, frg) in state.get('fragments'):
                self._fragments[k] = qclo.QcFragment(frg, qc_parent=self)
        self.charge = state.get('charge', 0)
        self._state = state.get('state', {})
        
    # pdfparam ---------------------------------------------------------
    def _get_pdfparam(self):
        '''
        pdfparamオブジェクトを返す
        '''
        pdfparam_path = os.path.abspath(os.path.join(self.work_dir,
                                                     self._pdfparam_filename))

        if 'pdfparam' not in self._cache:
            if os.path.exists(pdfparam_path):
                mpac_file = open(pdfparam_path, 'rb')
                mpac_data = msgpack.unpackb(mpac_file.read())
                mpac_data = bridge.Utils.to_unicode_dict(mpac_data)
                mpac_file.close()
                logger.debug('pdfparam({}) is loaded.'.format(pdfparam_path))
                self._cache['pdfparam'] = pdf.PdfParam(mpac_data)
            else:
                pdfsim = pdf.PdfSim()
                self._cache['pdfparam'] = pdf.get_default_pdfparam()
                logger.debug('use default pdfparam.')
        else:
            logger.debug('pdfparam is cached.')

        return self._cache['pdfparam']

    pdfparam = property(_get_pdfparam)
        

    # DB ---------------------------------------------------------------
    def _get_db_path(self):
        db_path = os.path.abspath(os.path.join(self.work_dir,
                                               self._db_filename))
        return db_path
    db_path = property(_get_db_path)

    def get_pdfarchive(self):
        '''
        pdfArchiveオブジェクトを返す
        '''
        pdfarc = pdf.PdfArchive(self.db_path)
        return pdfarc


    # ==================================================================
    # PROPERTIES
    # ==================================================================
    # work_dir ---------------------------------------------------------
    def _get_work_dir(self):
        return self._work_dir
    work_dir = property(_get_work_dir)

    # name -------------------------------------------------------------
    def _get_name(self):
        return self._name
    name = property(_get_name)

    # basisset ---------------------------------------------------------
    def _set_basisset(self, pdfparam):
        for fragment_name, fragment in self.fragments():
            fragment.set_basisset(pdfparam)

    # frame_molecule ---------------------------------------------------
    def _get_frame_molecule(self):
        '''
        モデリングされた分子構造をAtomGroupオブジェクトで返す
        '''
        if 'frame_molecule' not in self._cache:
            logger.info('create frame molecule coordinates.')
            frame_molecule = bridge.AtomGroup()
            for frg_name, frg in self._fragments.items():
                logger.info('fragment name={}: {} atoms'.format(frg_name,
                                                                      frg.get_number_of_all_atoms()))
                frame_molecule[frg_name] = frg.get_AtomGroup()
            self._cache['frame_molecule'] = frame_molecule
            logger.info('')
        return self._cache['frame_molecule']

    frame_molecule = property(_get_frame_molecule)

    # work dir ---------------------------------------------------------
    def _prepare_work_dir(self):
        '''
        カレントディレクトリ下に作業ディレクトリとなる
        nameのディレクトリを作成する。
        '''
        assert(len(self.name) > 0)
        self._work_dir = os.path.abspath(os.path.join(os.curdir, self.name))    
        if not os.path.exists(self.work_dir):
            logger.info('make work dir: {}'.format(self.work_dir))
            os.mkdir(self.work_dir)
        else:
            logger.debug('already exist: {}'.format(self.work_dir))

    def cd_work_dir(self, job_name=''):
        '''
        作業ディレクトリをオブジェクトのwork_dirに移動する
        '''
        logger.info('=' * 20)
        logger.info('>>>> {job_name}@{frame_name}'.format(job_name = job_name,
                                                                frame_name = self.name))
        logger.info('work dir: {work_dir}'.format(work_dir=self.work_dir))
        
        logger.info('=' * 20)
        self._prev_dir = os.path.abspath(os.curdir)
        os.chdir(self.work_dir)

    def restore_cwd(self):
        '''
        self.cd_work_dir() 以前のディレクトリに戻す
        '''
        os.chdir(self._prev_dir)
        logger.info('<<<< (prev_dir: {})\n'.format(self._prev_dir))

    def _check_path(self, path):
        if not os.path.exists(path):
            logger.warn('NOT FOUND: {}'.format(path))

    # charge -----------------------------------------------------------
    def _get_charge(self):
        return int(self._charge)

    def _set_charge(self, charge):
        self._charge = int(charge)

    charge = property(_get_charge, _set_charge)
            
    # num_of_AOs -------------------------------------------------------
    def get_number_of_AOs(self):
        '''
        return the number of atomic orbitals.
        '''
        num_of_AOs = 0
        for frg_name, frg in self.fragments():
            num_of_AOs += frg.get_number_of_AOs()
        return num_of_AOs

    # ==================================================================
    # STATE
    # ==================================================================
    # guess_density ----------------------------------------------------
    def _get_state_finished_guess_density(self):
        self._state.setdefault('is_finished_guess_density', False)
        return self._state['is_finished_guess_density']

    def _set_state_finished_guess_density(self, yn):
        self._state['is_finished_guess_density'] = bool(yn)
        
    is_finished_guess_density = property(_get_state_finished_guess_density,
                                         _set_state_finished_guess_density)
        
    # guess_QCLO -------------------------------------------------------
    def _get_state_finished_guess_QCLO(self):
        self._state.setdefault('is_finished_guess_QCLO', False)
        return self._state['is_finished_guess_QCLO']

    def _set_state_finished_guess_QCLO(self, yn):
        self._state['is_finished_guess_QCLO'] = bool(yn)
        
    is_finished_guess_QCLO = property(_get_state_finished_guess_QCLO,
                                      _set_state_finished_guess_QCLO)
        
    # pre-SCF ----------------------------------------------------------
    def _get_state_finished_prescf(self):
        self._state.setdefault('is_finished_prescf', False)
        return self._state['is_finished_prescf']

    def _set_state_finished_prescf(self, yn):
        self._state['is_finished_prescf'] = bool(yn)
        
    is_finished_prescf = property(_get_state_finished_prescf,
                                  _set_state_finished_prescf)
        
    # SCF --------------------------------------------------------------
    def _get_state_finished_scf(self):
        self._state.setdefault('is_finished_scf', False)
        return self._state['is_finished_scf']

    def _set_state_finished_scf(self, yn):
        self._state['is_finished_scf'] = bool(yn)
        
    is_finished_scf = property(_get_state_finished_scf,
                               _set_state_finished_scf)

    # Force ------------------------------------------------------------
    def _get_state_finished_force(self):
        self._state.setdefault('is_finished_force', False)
        return self._state['is_finished_force']

    def _set_state_finished_force(self, yn):
        self._state['is_finished_force'] = bool(yn)
        
    is_finished_force = property(_get_state_finished_force,
                                 _set_state_finished_force)
    
    
    # pick density matrix  ---------------------------------------------
    def _get_state_finished_pickup_density_matrix(self):
        self._state.setdefault('is_finished_pickup_density_matrix', False)
        return self._state['is_finished_pickup_density_matrix']

    def _set_state_finished_pickup_density_matrix(self, yn):
        self._state['is_finished_pickup_density_matrix'] = bool(yn)
        
    is_finished_pickup_density_matrix = property(_get_state_finished_pickup_density_matrix,
                                                 _set_state_finished_pickup_density_matrix)
    
    # LO ---------------------------------------------------------------
    def _get_state_finished_LO(self):
        self._state.setdefault('is_finished_LO', False)
        return self._state['is_finished_LO']

    def _set_state_finished_LO(self, yn):
        self._state['is_finished_LO'] = bool(yn)
        
    is_finished_LO = property(_get_state_finished_LO,
                               _set_state_finished_LO)
    
    # pickup LO --------------------------------------------------------
    def _get_state_finished_pickup_LO(self):
        self._state.setdefault('is_finished_pickup_LO', False)
        return self._state['is_finished_pickup_LO']

    def _set_state_finished_pickup_LO(self, yn):
        self._state['is_finished_pickup_LO'] = bool(yn)
        
    is_finished_pickup_LO = property(_get_state_finished_pickup_LO,
                                     _set_state_finished_pickup_LO)
    # ==================================================================
    # GUESS
    # ==================================================================
    # guess density ----------------------------------------------------
    def guess_density(self, run_type ='rks'):
        if self.is_finished_guess_density:
            logger.info('guess_density has been calced.')
            return

        self.cd_work_dir('guess_density')

        guess_density_matrix_path = 'guess.density.{}.mat'.format(run_type)

        # 既存のデータを消去する
        if os.path.exists(guess_density_matrix_path):
            os.remove(guess_density_matrix_path)
                
        pdfsim = pdf.PdfSim()
        pdfsim.setup()
        
        for frg_name, frg in self.fragments():
            logger.info('fragment name={}: {} atoms'.format(frg_name,
                                                                  frg.get_number_of_all_atoms()))
            if frg.qc_parent == None:
                logger.warn('guess_density(): qc_parent == None. frg_name={}'.format(frg_name))

            frg_guess_density_matrix_path = frg.prepare_guess_density_matrix(run_type)

            logger.debug('guess_density() [{}@{}] ext: {} from {}'.format(
                frg_name,
                frg.qc_parent.name,
                guess_density_matrix_path,
                frg_guess_density_matrix_path))
            if os.path.exists(frg_guess_density_matrix_path):
                pdf.run_pdf(['mat-ext', '-d',
                             guess_density_matrix_path,
                             frg_guess_density_matrix_path,
                             guess_density_matrix_path])
            else:
                logger.warn('not found: frg.guess.dens.mat={}'.format(frg_guess_density_matrix_path))

        self.pdfparam.guess = 'density_matrix'
        logger.info('initial guess (density matrix) created at {}'.format(guess_density_matrix_path))

        # check
        self._check_path(guess_density_matrix_path)
        
        self.is_finished_guess_density = True
        self._save()
            
        self.restore_cwd()

        
    def guess_QCLO(self, run_type='rks', isCalcOrthogonalize = False):
        if self.is_finished_guess_QCLO:
            logger.info('guess_density has been calced.')
            return

        if self.is_finished_prescf != True:
            self.calc_preSCF()

        self.cd_work_dir('guess_QCLO')

        guess_QCLO_matrix_path = 'guess.QCLO.{}.mat'.format(run_type)
        if os.path.exists(guess_QCLO_matrix_path):
            os.remove(guess_QCLO_matrix_path)
        
        num_of_AOs = 0
        for frg_name, frg in self.fragments():
            logger.info('guess QCLO: frg_name={}, parent={}'.format(frg_name, frg.qc_parent.name))
            frg_QCLO_matrix_path = frg.prepare_guess_QCLO_matrix(run_type, self)
            if os.path.exists(frg_QCLO_matrix_path):
                pdf.run_pdf(['mat-ext', '-c',
                             guess_QCLO_matrix_path,
                             frg_QCLO_matrix_path,
                             guess_QCLO_matrix_path])
            else:
                logger.warn('The QCLO of the subgroup, {}, was not created.'.format(frg_name))

        # orthogonalize
        guess_path = 'guess.lcao.{}.mat'.format(run_type)
        if isCalcOrthogonalize:
            logger.info('orthogonalize')
            Xinv_path = self.pdfparam.get_Xinv_mat_path()

            self._check_path(guess_QCLO_matrix_path)
            pdf.run_pdf(['mat-mul', '-v',
                         Xinv_path,
                         guess_QCLO_matrix_path,
                         guess_path])
        else:
            shutil.copy(guess_QCLO_matrix_path, guess_path)

        self.pdfparam.guess = 'lcao'
        logger.info('guess LCAO matrix created: {}'.format(guess_path))

        # check
        self._check_path(guess_QCLO_matrix_path)

        
        self.is_finished_guess_QCLO = True
        self._save()
        self.restore_cwd()

        # create occ file
        self._create_occupation_file(run_type)

    def _create_occupation_file(self, run_type='rks'):
        self.cd_work_dir('create occ')

        occ_level = -1
        electrons_per_orb = 0.0
        run_type = run_type.upper()
        if run_type == 'RKS':
            occ_level = int((self.pdfparam.num_of_electrons / 2.0))
            electrons_per_orb = 2.0
        else:
            logger.critical('NOT supported. run_type={}'.format(run_type))

        num_of_MOs = self.pdfparam.num_of_MOs
        occ_vtr = pdf.Vector(num_of_MOs)
        for i in range(occ_level):
            occ_vtr.set(i, electrons_per_orb)

        occ_vtr_path = 'guess.occ.{}.vtr'.format(run_type.lower())
        occ_vtr.save(occ_vtr_path)
        self._check_path(occ_vtr_path)

        self._save()
        self.restore_cwd()
        
    # ==================================================================
    # CALC
    # ==================================================================
    # ------------------------------------------------------------------
    def calc_preSCF(self, dry_run=False):
        '''
        '''
        if self.is_finished_prescf:
            logger.info('preSCF has been calced.')
            return
        
        self.cd_work_dir('calc preSCF')

        pdfsim = pdf.PdfSim()
        for frg_name, frg in self.fragments():
            frg.set_basisset(self.pdfparam)
        self.pdfparam.molecule = self.frame_molecule

        # num_of_electrons
        num_of_electrons = self.pdfparam.num_of_electrons # calc from the molecule data
        logger.info('the number of electrons = {}'.format(num_of_electrons))
        if self.charge != 0:
            logger.info('specify the charge => {}'.format(self.charge))
            num_of_electrons -= self.charge # 電子(-)数と電荷(+)の正負が逆なことに注意
            self.pdfparam.num_of_electrons = num_of_electrons
            logger.info('update the number of electrons => {}'.format(self.pdfparam.num_of_electrons))
        
        self.pdfparam.step_control = 'integral'
        pdfsim.sp(self.pdfparam,
                  workdir = self.work_dir,
                  db_path = self.db_path,
                  dry_run = dry_run)
        
        self._cache.pop('pdfparam')
        self.is_finished_prescf = True
        self._save()

        self.restore_cwd()
        
    # sp ---------------------------------------------------------------
    def calc_sp(self, dry_run=False):
        '''
        calculate single point energy
        '''
        if self.is_finished_scf:
            logger.info('SP has been calced.')
            return
            
        if self.is_finished_prescf != True:
            self.calc_preSCF(dry_run)

        self.cd_work_dir('calc SP')

        pdfsim = pdf.PdfSim()
        for frg_name, frg in self.fragments():
            frg.set_basisset(self.pdfparam)
        self.pdfparam.molecule = self.frame_molecule

        # num_of_electrons
        num_of_electrons = self.pdfparam.num_of_electrons # calc from the molecule data
        logger.info('the number of electrons = {}'.format(num_of_electrons))
        if self.charge != 0:
            logger.info('specify the charge => {}'.format(self.charge))
            num_of_electrons -= self.charge # 電子(-)数と電荷(+)の正負が逆なことに注意
            self.pdfparam.num_of_electrons = num_of_electrons
            logger.info('update the number of electrons => {}'.format(self.pdfparam.num_of_electrons))
        
        #self.output_xyz("{}/model.xyz".format(self.name))

        self.pdfparam.step_control = 'guess scf'
        pdfsim.sp(self.pdfparam,
                  workdir = self.work_dir,
                  db_path = self.db_path,
                  dry_run = dry_run)

        self._cache.pop('pdfparam')
        self.is_finished_scf = True
        self._grouping_fragments()
        self._switch_fragments()
        self._save()
            
        self.restore_cwd()

    # force ------------------------------------------------------------
    def calc_force(self, dry_run=False):
        '''
        calculate force (energy gradient)

        absolute: force -> gradient
        '''
        if self.is_finished_force:
            logger.info('force has been calced.')
            return
            
        if self.is_finished_scf != True:
            self.calc_sp(dry_run)

        self.cd_work_dir('calc force')

        pdfsim = pdf.PdfSim()
        for frg_name, frg in self.fragments():
            frg.set_basisset(self.pdfparam)
        self.pdfparam.molecule = self.frame_molecule

        # num_of_electrons
        num_of_electrons = self.pdfparam.num_of_electrons # calc from the molecule data
        logger.info('the number of electrons = {}'.format(num_of_electrons))
        if self.charge != 0:
            logger.info('specify the charge => {}'.format(self.charge))
            num_of_electrons -= self.charge # 電子(-)数と電荷(+)の正負が逆なことに注意
            self.pdfparam.num_of_electrons = num_of_electrons
            logger.info('update the number of electrons => {}'.format(self.pdfparam.num_of_electrons))

        self.pdfparam.step_control = 'force'
        pdfsim.sp(self.pdfparam,
                  workdir = self.work_dir,
                  db_path = self.db_path,
                  dry_run = dry_run)

        self._cache.pop('pdfparam')
        self.is_finished_force = True
        self._save()
            
        self.restore_cwd()

    # summary ------------------------------------------------------------------
    def summary(self, dry_run=False, format_str=None, filepath=None):
        '''
        Format: 
            {NUM_OF_ATOMS}: number of atoms
            {NUM_OF_AO}:    number of AOs
            {NUM_OF_MO}:    number of MOs
            {METHOD}:       method
            {IS_CONVERGED}: Whether the SCF is converged or not
            {ITERATION}:    iteration
            {TOTAL_ENERGY}: total energy
            {GRADIENT_RMS}: gradient RMS
        '''
        if self.is_finished_scf != True:
            self.calc_sp(dry_run)

        self.cd_work_dir('summary')

        values = {}
        pdfarc = self.get_pdfarchive()
        values['NUM_OF_ATOMS'] = pdfarc.num_of_atoms
        values['NUM_OF_AO'] = pdfarc.num_of_AOs
        values['NUM_OF_MO'] = pdfarc.num_of_MOs
        values['METHOD'] = pdfarc.method
        values['IS_CONVERGED'] = pdfarc.scf_converged
        itr = pdfarc.iterations
        values['ITERATION'] = itr
        values['TOTAL_ENERGY'] = pdfarc.get_total_energy(itr)
        values['GRADIENT_RMS'] = pdfarc.get_gradient_rms()

        if format_str == None:
            format_str = 'total energy: {TOTAL_ENERGY} at {ITERATION}'
        output = format_str.format(**values)

        if output[-1] != "\n":
            output += "\n"
        
        logger.info(output)
        if filepath != None:
            with open(filepath, 'a') as f:
                f.write(output)

        self.restore_cwd()
        return output

    
    def get_gradient(self):
        '''
        '''
        self.cd_work_dir('get_gradient')

        pdfarc = self.get_pdfarchive()
        num_of_atoms = pdfarc.num_of_atoms
        grad =[ [] * num_of_atoms]
        for atom_index in range(num_of_atoms):
            grad[atom_index] = pdfarc.get_force(atom_index)
        
        self.restore_cwd()
    
    # ------------------------------------------------------------------
    

    # ==================================================================
    # PICKUP
    # ==================================================================
    # pickup density matrix --------------------------------------------
    def pickup_density_matrix(self, runtype ='rks'):
        '''
        密度行列を各フラグメントに割り当てる
        '''
        if self.is_finished_pickup_density_matrix:
            logger.info('pickup density matrix has done.')
            return
        
        self.cd_work_dir('pickup density matrix')
        
        dens_mat_path = self.pdfparam.get_density_matrix_path(runtype=runtype)
        logger.info('ref. density matrix: {}'.format(dens_mat_path))
        
        global_dim = 0
        for frg_name, frg in self.fragments():
            dim = frg.get_number_of_AOs()
            if dim > 0:
                frg_dens_mat_path = 'Ppq.{}.{}.mat'.format(runtype, frg_name)
                logger.info(
                    'select [{start}:{end}] for [{frame}][{fragment}]'.format(
                        frame=frg_name,
                        fragment=self.name,
                        start=global_dim,
                        end=global_dim +dim -1))

                # フラグメント対応部分を切り出す
                pdf.run_pdf(['mat-select',
                             '-t', global_dim,
                             '-l', global_dim,
                             '-b', global_dim +dim -1,
                             '-r', global_dim +dim -1,
                             dens_mat_path,
                             frg_dens_mat_path])

                # select された行列を対称行列に変換
                pdf.run_pdf(['mat-symmetrize',
                             frg_dens_mat_path,
                             frg_dens_mat_path])

                # 対称行列のパスをフラグメントに登録
                frg.set_density_matrix(frg_dens_mat_path)

                logger.info(
                    'density matrix for [{frame}][{fragment}] was saved as {path}'.format(
                        fragment=frg_name,
                        frame=self.name,
                        path=frg_dens_mat_path))
                global_dim += dim

        logger.is_finished_pickup_density_matrix = True
        self._save()
        self.restore_cwd()
            
    # ------------------------------------------------------------------
    def calc_lo(self, run_type, dry_run=False):
        if self.is_finished_LO:
            logger.info('LO has done.')
            return

        if self.is_finished_scf != True:
            self.calc_sp(dry_run=dry_run)
        
        self.cd_work_dir('calc lo')

        logger.info('start lo calculation.')
        pdf.run_pdf('lo')

        self.is_finished_LO = True
        self._save()
        self.restore_cwd()

    # ------------------------------------------------------------------
    def pickup_QCLO(self, run_type='rks'):
        if self.is_finished_pickup_LO:
            logger.info('pickup LO has been finished.')
            return

        self.calc_lo(run_type)
            
        self.cd_work_dir('pickup lo')
        
        # debug
        pdfarc = self.get_pdfarchive()
        num_of_AOs = pdfarc.num_of_AOs
        num_of_MOs = pdfarc.num_of_MOs
        HOMO_level = pdfarc.get_HOMO_level('rks') # option base 0
        logger.info('num of AOs: {}'.format(num_of_AOs))
        logger.info('num of MOs: {}'.format(num_of_MOs))
        logger.info('HOMO level: {}'.format(HOMO_level +1)) 

        logger.info('fragment information:')
        for frg_name, frg in self.fragments():
            frg_AOs = frg.get_number_of_AOs()
            logger.info('fragment name:[{}] AOs={}'.format(frg_name, frg_AOs))
        logger.info('')
            
        # calc S*C
        if 'pdfparam' in self._cache:
            self._cache.pop('pdfparam')
        lo_satisfied = self.pdfparam.lo_satisfied
        if lo_satisfied != True:
            logger.warn('lo_satisfied: {}'.format(lo_satisfied))
        lo_iterations = self.pdfparam.lo_num_of_iterations
        logger.info('lo iterations: {}'.format(lo_iterations))

        logger.info('calc S*C')
        CSC_path = 'CSC.mat'
        Clo_path = self.pdfparam.get_clomat_path()
        pdf.run_pdf(['component',
                     '-v',
                     '-S', 'CSC.mat',
                     '-c', Clo_path])

        # load CSC
        CSC = pdf.Matrix()
        CSC.load(CSC_path)

        logger.info('make AO v.s. fragment table')
        AO_frg_tbl = self._get_AO_fragment_table(num_of_AOs)
        
        # pickup
        logger.info('assign fragment: start: HOMO={}'.format(HOMO_level))
        MO_fragment_assigned = {}        
        for mo in range(HOMO_level +1):
            frg_name = self._define_lo_fragment(mo, num_of_AOs, AO_frg_tbl, CSC)
            logger.info('{}-th MO -> fragment: [{}]'.format(mo, frg_name))
            MO_fragment_assigned.setdefault(frg_name, [])
            MO_fragment_assigned[frg_name].append(mo)
        logger.info('assign fragment: end')

        # assign report
        logger.info('==== assign report ====')
        for k, MOs in MO_fragment_assigned.items():
            logger.info('fragment[{}] # of MOs={}'.format(k, len(MOs)))
        logger.info('')
            
        # フラグメントのC_LOを作成する
        logger.info('create C_LO: start')
        Clo = pdf.Matrix()
        Clo.load(Clo_path)
        assert(num_of_AOs == Clo.rows)
        for frg_name, frg in self.fragments():
            frg_cols = len(MO_fragment_assigned.get(frg_name, []))
            logger.info('frg {}: col={}'.format(frg_name, frg_cols))
            if frg_cols == 0:
                logger.warn('fragment "{}" has no colomns.'.format(frg_name))
                # continue

            Clo_frg = pdf.Matrix(num_of_AOs, frg_cols)
            if frg_name in MO_fragment_assigned:
                for col, ref_col in enumerate(MO_fragment_assigned[frg_name]):
                    for row in range(num_of_AOs):
                        v = Clo.get(row, ref_col)
                        Clo_frg.set(row, col, v)

            Clo_path = 'Clo_{}.mat'.format(frg_name)
            logger.info('fragment C_LO save: {}'.format(Clo_path))
            Clo_frg.save(Clo_path)
            frg.set_LO_matrix(Clo_path, run_type)
        logger.info('create C_LO: end')
            
        # trans C_LO
        self._trans_LO()

        # finish
        self.is_finished_pickup_LO = True
        self._save()
        self.restore_cwd()

    def _get_AO_fragment_table(self, num_of_AOs):
        '''
        AO v.s. fragment_name の辞書を返す
        '''
        frg_table = [ None for x in range(num_of_AOs) ]

        AO_index = 0
        for frg_name, frg in self.fragments():
            frg_num_of_AOs = frg.get_number_of_AOs()
            for i in range(AO_index, AO_index + frg_num_of_AOs):
                frg_table[i] = frg_name
            AO_index += frg_num_of_AOs

        return frg_table

    def _define_lo_fragment(self, mo, num_of_AOs, AO_frg_tbl, CSC):
        judge = {}
        total = 0.0
        for ao in range(num_of_AOs):
            frg_name = AO_frg_tbl[ao]
            v = math.fabs(CSC.get(ao, mo))
            total += v
            judge.setdefault(frg_name, 0.0)
            judge[frg_name] += v
            
        for frg_name in judge.keys():
            judge[frg_name] /= total

        ranked_judge = sorted(judge.items(), key=lambda x:x[1], reverse=True)

        for rank, (k, v) in enumerate(ranked_judge):
            logger.info('[{rank}] name:{name}, score:{score:.3f}'.format(
                rank=rank +1,
                name=k,
                score=v))

        high_score = ranked_judge[0][1]
        if high_score < 0.5:
            logger.warn('1st score is too small: {}'.format(high_score))
            
        return ranked_judge[0][0]

    def _trans_LO(self):
        logger.info('trans LO at {}'.format(os.getcwd()))
        run_type = 'rks'
        F_path = self.pdfparam.get_Fmat_path(run_type)
        logger.info('F matrix: {}'.format(F_path))
        for frg_name, frg in self.fragments():
            logger.info('frg {}: AOs={}'.format(frg_name, frg.get_number_of_AOs()))
            if frg.get_number_of_AOs() == 0:
                continue
            
            Clo_path = frg.get_LO_matrix_path(run_type)
            #if Clo_path is None:
            #    continue

            # calc (C_LO)dagger * F * C_LO => F'
            F_Clo_path = 'F_Clo.{}.mat'.format(frg_name)
            pdf.run_pdf(['mat-mul', '-v',
                         F_path,
                         Clo_path,
                         F_Clo_path])

            Clo_dagger_path = 'Clo_dagger.{}.mat'.format(frg_name)
            pdf.run_pdf(['mat-transpose', '-v',
                         Clo_path,
                         Clo_dagger_path])

            F_prime_path = 'Fprime.{}.mat'.format(frg_name)
            pdf.run_pdf(['mat-mul', '-v',
                         Clo_dagger_path,
                         F_Clo_path,
                         F_prime_path])

            pdf.run_pdf(['mat-symmetrize',
                         F_prime_path,
                         F_prime_path])
            
            # diagonal F'
            eigval_path = 'QCLO_eigval.{}.vtr'.format(frg_name)
            Cprime_path = 'Cprime.{}.mat'.format(frg_name)
            logger.info("diagonal F'")
            pdf.run_pdf(['mat-diagonal', '-v',
                         '-l', eigval_path,
                         '-x', Cprime_path,
                         F_prime_path])

            # AO基底に変換
            C_QCLO_path = 'C_QCLO.{}.mat'.format(frg_name)
            pdf.run_pdf(['mat-mul', '-v',
                         Clo_path,
                         Cprime_path,
                         C_QCLO_path])
            frg.set_QCLO_matrix(C_QCLO_path)

            logger.info('C_QCLO saved: {}'.format(C_QCLO_path))
        
        
    #  =================================================================
    #  for fragments
    #  =================================================================
    def fragments(self):
        '''
        フラグメントの名前とオブジェクトを返すイテレータ
        '''
        for k in self._fragments.keys():
            yield(k, self._fragments[k])

    def has_fragment(self, fragment_name):
        '''
        フラグメントを持っていればTrueを返す
        '''
        fragment_name = bridge.Utils.to_unicode(fragment_name)

        return fragment_name in self._fragments.keys()
            
    # operator[] -------------------------------------------------------
    def __getitem__(self, fragment_name):
        '''
        出力用[]演算子
        '''
        fragment_name = bridge.Utils.to_unicode(fragment_name)
        return self._fragments.get(fragment_name, None)

    def __setitem__(self, fragment_name, obj):
        '''
        入力用[]演算子

        計算前であれば代入可能(つまりモデリング中)であるが、
        計算後は代入できない
        '''
        if self.is_finished_scf:
            logger.warn('operator[] called after simulation.')
            return

        if 'frame_molecule' in self._cache:
            self._cache.pop('frame_molecule')

        fragment_name = bridge.Utils.to_unicode(fragment_name)

        if isinstance(obj, qclo.QcFragment):
            fragment = qclo.QcFragment(obj)
            fragment.name = fragment_name
            if fragment.qc_parent == None:
                fragment.qc_parent = self
            logger.info('add fragment: name={}'.format(fragment_name))
            self._fragments[fragment_name] = fragment
        elif isinstance(obj, qclo.QcFrame):
            logger.info('add frame molecule: for {}'.format(fragment_name))
            fragment = qclo.QcFragment()
            fragment.name = fragment_name
            for k, f in obj.fragments():
                if not f.margin:
                    logger.warn('add fragment: k={} for {}'.format(k, fragment_name))
                    fragment.set_group(k, f)
                else:
                    logger.warn('pass fragment: k={} is margin'.format(k))
            if fragment.qc_parent == None:
                fragment.qc_parent = self
            self._fragments[fragment_name] = fragment
        else:
            raise

    # rearangement -----------------------------------------------------
    def _switch_fragments(self):
        '''
        fragmentsを入力用から出力用に切り替える

        処理内容:
        - 各fragmentの親を自分(self)にする
        '''
        logger.info('switch fragment.')
        output_fragments = OrderedDict()
        for frg_name, frg in self.fragments():
            logger.info('fragment_name: {}'.format(frg_name))
            new_frg = qclo.QcFragment(frg, qc_parent=self)
            assert(new_frg.qc_parent.name == self.name)
            output_fragments[frg_name] = new_frg
        self._fragments = output_fragments

        #logger.info('merge subgroups')
        #for key, frg in self.fragments():
        #    frg.merge_subgroups()
        
        logger.info('---> switch ')
        for frg_name, frg in self.fragments():
            logger.info('{}: parent={}'.format(frg_name,
                                                     frg.qc_parent.name))
        logger.info('<---')

    def _grouping_fragments(self):
        logger.info('grouping fragments.')
        for frg_name, frg in self.fragments():
            frg.grouping_subfragments()

    # ==================================================================
    # coordinates
    # ==================================================================
    # outout XYZ -------------------------------------------------------
    def output_xyz(self, file_path):
        xyz = bridge.Xyz(self.frame_molecule)
        xyz.save(file_path)


    # ==================================================================
    # orbital table
    # ==================================================================
    def get_orbital_info(self):
        '''
        AOに対するQcOrbitalDataリストを返す
        '''
        orbinfo = []
        for k, frg in self.fragments():
            orbinfo.extend(frg.get_orbital_info())
        return orbinfo
        
    
    # ==================================================================
    # operators
    # ==================================================================
    # operator == ------------------------------------------------------
    def __eq__(self, rhs):
        if rhs == None:
            return False
        return (self.name == rhs.name)
        
    def __ne__(self, rhs):
        return not self.__eq__(rhs)

    # operator str -----------------------------------------------------
    def __str__(self):
        answer = ""
        answer = 'frame name={}\n'.format(self.name)
        for key, fragment in self.fragments():
            answer += '>> fragment: {}\n'.format(key)
            answer += str(fragment)
            answer += '\n'
        return answer

