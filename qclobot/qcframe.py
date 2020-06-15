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

from .qcfragment import QcFragment
import proteindf_tools as pdf
import proteindf_bridge as bridge
import shutil
from collections import OrderedDict
import math
import os
import logging
logger = logging.getLogger(__name__)

try:
    import msgpack
except:
    import msgpack_pure as msgpack


class QcFrame(object):
    _pdfparam_filename = 'pdfparam.mpac'
    _db_filename = 'pdfresults.h5'
    TOO_SMALL = 1.0E-5

    # ------------------------------------------------------------------
    def __init__(self, name, *args, **kwargs):
        '''
        create QcFrame object.

        name: name of the frame molecule
        '''
        # mandatory parameter
        self._name = name
        self._fragments = OrderedDict()
        self._charge = 0
        self._state = {}  # 状態保持
        self._cmds = self._get_default_cmds()  # 外部コマンド

        self._initialize()

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
        self._cmds = copy.copy(rhs._cmds)

    # def __del__(self):
    #    self._save()

    def _initialize(self, *args, **kwargs):
        pass

    def _get_default_cmds(self):
        answer = {}
        answer['mat-extend'] = 'mat-extend'
        answer['mat-mul'] = 'mat-mul'
        answer['mat-select'] = 'mat-select'
        answer['mat-symmetrize'] = 'mat-symmetrize'
        answer['mat-transpose'] = 'mat-transpose'
        answer['mat-diagonal'] = 'mat-diagonal'
        answer['archive'] = 'archive-h5'

        return answer

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

    def save(self):
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
        state['cmds'] = self._cmds
        return state

    def __setstate__(self, state):
        assert(isinstance(state, dict))
        self._name = state.get('name')
        self._fragments = OrderedDict()
        if 'fragments' in state:
            for (k, frg) in state.get('fragments'):
                self._fragments[k] = QcFragment(frg, parent=self)
        self.charge = state.get('charge', 0)
        self._state = state.get('state', {})
        self._cmds = state.get('cmds', self._get_default_cmds())

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
    def set_db_filename(self, filename):
        self._db_filename = str(filename)

    def _get_db_path(self):
        db_path = os.path.abspath(os.path.join(self.work_dir,
                                               self._db_filename))
        return db_path
    db_path = property(_get_db_path)

    def get_pdfarchive(self):
        '''
        pdfArchiveオブジェクトを返す
        '''
        logger.debug("get_pdfarchive db_path={}".format(self.db_path))
        pdfarc = None
        if self._cmds.get('archive', None) == 'archive':
            pdfarc = pdf.PdfArchive(self.db_path)
        else:
            pdfarc = pdf.PdfParam_H5(self.db_path)
        return pdfarc

    # ==================================================================
    # PROPERTIES
    # ==================================================================
    # command alias ----------------------------------------------------

    def set_command_alias(self, cmd_alias_dict):
        for k, v in cmd_alias_dict.items():
            logger.debug("command update: {} -> {}".format(k, v))
            self._cmds[k] = v

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
                logger.info('fragment name={name}: atoms={atoms}, elec={elec}, charge={charge}'.format(
                    name=frg_name,
                    atoms=frg.get_number_of_all_atoms(),
                    elec=frg.sum_of_atomic_number(),
                    charge=frg.get_AtomGroup().charge))
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
        # assert(len(self.name) > 0)
        if len(self.name) == 0:
            logger.critical("frame name is not defined.")
            raise
        self._work_dir = os.path.abspath(os.path.join(os.curdir, self.name))
        if not os.path.exists(self.work_dir):
            logger.info("{header} make work dir: {path}".format(
                header=self.header,
                path=self.work_dir))
            os.mkdir(self.work_dir)
        else:
            logger.debug("{header} already exist: {path}".format(
                header=self.header,
                path=self.work_dir))

    def cd_work_dir(self, job_name=''):
        '''
        作業ディレクトリをオブジェクトのwork_dirに移動する
        '''
        logger.info('=' * 20)
        logger.info("{header} > {job_name}@{frame_name}".format(
            header=self.header,
            job_name=job_name,
            frame_name=self.name))
        logger.debug("{header} work dir: {work_dir}".format(
            header=self.header,
            work_dir=self.work_dir))

        logger.info('=' * 20)
        self._prev_dir = os.path.abspath(os.curdir)
        os.chdir(self.work_dir)

    def restore_cwd(self):
        '''
        self.cd_work_dir() 以前のディレクトリに戻す
        '''
        os.chdir(self._prev_dir)
        logger.debug("{header} < (prev_dir: {path})".format(
            header=self.header,
            path=self._prev_dir))

    def _check_path(self, path):
        if not os.path.exists(path):
            logger.warn("{header} NOT FOUND: {path}".format(
                header=self.header, path=path))

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

    def guess_density(self, run_type='rks', force=False):
        if ((self.is_finished_guess_density == True) and
                (force == False)):
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
            if frg.parent == None:
                logger.warn(
                    'guess_density(): parent == None. frg_name={}'.format(frg_name))

            frg.set_command_alias(self._cmds)
            frg_guess_density_matrix_path = frg.prepare_guess_density_matrix(
                run_type)

            logger.debug('guess_density() [{}@{}] ext: {} from {}'.format(
                frg_name,
                frg.parent.name,
                guess_density_matrix_path,
                frg_guess_density_matrix_path))
            if os.path.exists(frg_guess_density_matrix_path):
                pdf.run_pdf([self._cmds['mat-extend'], '-d',
                             guess_density_matrix_path,
                             frg_guess_density_matrix_path,
                             guess_density_matrix_path])
            else:
                logger.warn('not found: frg.guess.dens.mat={}'.format(
                    frg_guess_density_matrix_path))

        self.pdfparam.guess = 'density_matrix'
        logger.info('initial guess (density matrix) created at {}'.format(
            guess_density_matrix_path))

        # check
        self._check_path(guess_density_matrix_path)

        self.is_finished_guess_density = True
        self.save()

        self.restore_cwd()

    def guess_QCLO(self, run_type='rks',
                   force=False,
                   isCalcOrthogonalize=False):
        """create guess by using QCLO method
        """
        if ((self.is_finished_guess_QCLO == True) and
                (force == False)):
            logger.info('guess_density has been calced.')
            return

        self.cd_work_dir('guess_QCLO')

        guess_QCLO_matrix_path = 'guess.QCLO.{}.mat'.format(run_type)
        if os.path.exists(guess_QCLO_matrix_path):
            os.remove(guess_QCLO_matrix_path)

        num_of_AOs = 0
        for frg_name, frg in self.fragments():
            logger.info('guess QCLO: frg_name={}, parent={}'.format(
                frg_name, frg.parent.name))

            frg.set_command_alias(self._cmds)
            frg_QCLO_matrix_path = frg.prepare_guess_QCLO_matrix(
                run_type, self, force=force)
            if os.path.exists(frg_QCLO_matrix_path):
                pdf.run_pdf([self._cmds['mat-extend'], '-c',
                             guess_QCLO_matrix_path,
                             frg_QCLO_matrix_path,
                             guess_QCLO_matrix_path])
            else:
                logger.warn(
                    'The QCLO of the subgroup, {}, was not created.'.format(frg_name))

        # orthogonalize
        guess_path = 'guess.lcao.{}.mat'.format(run_type)
        if isCalcOrthogonalize:
            if self.is_finished_prescf != True:
                self.calc_preSCF()

            logger.info('orthogonalize')
            Xinv_path = self.pdfparam.get_Xinv_mat_path()

            self._check_path(guess_QCLO_matrix_path)
            pdf.run_pdf([self._cmds['mat-mul'], '-v',
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
        self.save()
        self.restore_cwd()

        # create occ file
        self._create_occupation_file(run_type)

    def _create_occupation_file(self, run_type='rks'):
        self.cd_work_dir('create occ')

        self._setup_pdf()

        occ_level = -1
        electrons_per_orb = 0.0
        run_type = run_type.upper()
        if run_type == 'RKS':
            occ_level = int((self.pdfparam.num_of_electrons / 2.0))
            electrons_per_orb = 2.0
        else:
            logger.critical("{header} NOT supported. run_type={run_type}".format(
                header=self.header, run_type=run_type))

        # num_of_MOs = self.pdfparam.num_of_MOs
        # occ_vtr = pdf.Vector(num_of_MOs)
        occ_vtr = pdf.Vector(occ_level)
        for i in range(occ_level):
            occ_vtr.set(i, electrons_per_orb)

        occ_vtr_path = "guess.occ.{}.vtr".format(run_type.lower())
        occ_vtr.save(occ_vtr_path)
        self._check_path(occ_vtr_path)

        self.save()
        self.restore_cwd()

    # ==================================================================
    # CALC
    # ==================================================================

    def _setup_pdf(self):
        logger.info("{header} setup ProteinDF condition".format(
            header=self.header))

        for frg_name, frg in self.fragments():
            frg.set_basisset(self.pdfparam)
        self.pdfparam.molecule = self.frame_molecule

        # num_of_electrons
        # calc from the molecule data
        num_of_electrons = self.frame_molecule.sum_of_atomic_number()
        logger.info("{header} the number of electrons = {elec}".format(
            header=self.header, elec=num_of_electrons))
        if self.charge != 0:
            logger.info('specify the charge => {}'.format(self.charge))
            num_of_electrons -= self.charge  # 電子(-)数と電荷(+)の正負が逆なことに注意
            self.pdfparam.num_of_electrons = num_of_electrons
            logger.info("{header} update the number of electrons => {elec}".format(
                header=self.header,
                elec=self.pdfparam.num_of_electrons))
        if self.pdfparam.num_of_electrons % 2 != 0:
            logger.warning(
                "{header} the number of electrons is not even.".format(header=self.header))

    # ------------------------------------------------------------------
    def calc_preSCF(self, dry_run=False):
        '''
        '''
        if self.is_finished_prescf:
            logger.info('preSCF has been calced.')
            return

        self.cd_work_dir('calc preSCF')
        self.check_bump_of_atoms()

        self._setup_pdf()
        self.pdfparam.step_control = 'integral'
        self.save()

        pdfsim = pdf.PdfSim()
        pdfsim.sp(self.pdfparam,
                  workdir=self.work_dir,
                  db_path=self.db_path,
                  dry_run=dry_run,
                  cmd_archive=self._cmds['archive'])

        self._cache.pop('pdfparam')
        self.is_finished_prescf = True
        self.save()

        self.restore_cwd()

    # sp ---------------------------------------------------------------
    def calc_sp(self, dry_run=False):
        '''
        calculate single point energy
        '''
        if self.is_finished_scf:
            logger.info('SP has been calced.')
            self._grouping_fragments()
            self._switch_fragments()
            return

        if self.is_finished_prescf != True:
            self.calc_preSCF(dry_run)

        self.cd_work_dir('calc SP')
        self.check_bump_of_atoms()

        self._setup_pdf()
        # self.output_xyz("{}/model.xyz".format(self.name))
        self.pdfparam.step_control = 'guess scf'
        self.save()

        pdfsim = pdf.PdfSim()
        pdfsim.sp(self.pdfparam,
                  workdir=self.work_dir,
                  db_path=self.db_path,
                  dry_run=dry_run,
                  cmd_archive=self._cmds['archive'])

        self._cache.pop('pdfparam')
        self.is_finished_scf = True

        self._grouping_fragments()
        self._switch_fragments()
        self.save()

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

        self._setup_pdf()
        self.pdfparam.step_control = 'force'
        self.save()

        pdfsim = pdf.PdfSim()
        # for frg_name, frg in self.fragments():
        #     frg.set_basisset(self.pdfparam)
        # self.pdfparam.molecule = self.frame_molecule
        #
        # # num_of_electrons
        # num_of_electrons = self.pdfparam.num_of_electrons # calc from the molecule data
        # logger.info('the number of electrons = {}'.format(num_of_electrons))
        # if self.charge != 0:
        #     logger.info('specify the charge => {}'.format(self.charge))
        #     num_of_electrons -= self.charge # 電子(-)数と電荷(+)の正負が逆なことに注意
        #     self.pdfparam.num_of_electrons = num_of_electrons
        #     logger.info('update the number of electrons => {}'.format(self.pdfparam.num_of_electrons))
        pdfsim.sp(self.pdfparam,
                  workdir=self.work_dir,
                  db_path=self.db_path,
                  dry_run=dry_run,
                  cmd_archive=self._cmds['archive'])

        self._cache.pop('pdfparam')
        self.is_finished_force = True
        self.save()

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
        grad = [[] * num_of_atoms]
        for atom_index in range(num_of_atoms):
            grad[atom_index] = pdfarc.get_force(atom_index)

        self.restore_cwd()

    # pop --------------------------------------------------------------
    def pop(self, dry_run=False, iteration=-1):
        '''
        '''
        if self.is_finished_scf != True:
            self.calc_sp(dry_run)

        if iteration == -1:
            iteration = self.pdfparam.iterations

        self._calc_pop(iteration=iteration)
        pop_vtr = self.get_pop(iteration)

        self.save()
        self.restore_cwd()

        return pop_vtr

    def _calc_pop(self, iteration=-1, dry_run=False):
        """
        """
        if iteration == -1:
            iteration = self.pdfparam.iterations
        self.cd_work_dir('calc pop: iteration={}'.format(iteration))

        pdfsim = pdf.PdfSim()
        pdfsim.pop(iteration=iteration,
                   dry_run=dry_run)

        self.restore_cwd()

    def get_pop(self, iteration=-1):
        """
        """
        if iteration == -1:
            iteration = self.pdfparam.iterations
        self.cd_work_dir('get pop: iteration={}'.format(iteration))

        run_type = "rks"
        pop_path = self.pdfparam.get_pop_mulliken_path(
            run_type, iteration=iteration)
        pop_vtr = pdf.Vector()
        pop_vtr.load(pop_path)

        self.restore_cwd()

        return pop_vtr

    # ------------------------------------------------------------------

    # ==================================================================
    # PICKUP
    # ==================================================================
    # pickup density matrix --------------------------------------------

    def pickup_density_matrix(self, runtype='rks'):
        '''
        密度行列を各フラグメントに割り当てる
        '''
        if self.is_finished_pickup_density_matrix:
            logger.info("{header} pickup density matrix has done.".format(
                header=self.header))
            return

        self.cd_work_dir('pickup density matrix')

        # post-SCF
        self._grouping_fragments()
        self._switch_fragments()

        dens_mat_path = self.pdfparam.get_density_matrix_path(runtype=runtype)
        logger.info("{header} reference density matrix: {path}".format(
            header=self.header,
            path=dens_mat_path))

        global_dim = 0
        for frg_name, frg in self.fragments():
            dim = frg.get_number_of_AOs()
            if dim > 0:
                frg_dens_mat_path = 'Ppq.{}.{}.mat'.format(runtype, frg_name)
                logger.info("{header} select [{start}:{end}] for {fragment}".format(
                    header=self.header,
                    fragment=frg_name,
                    start=global_dim,
                    end=global_dim + dim - 1))

                # フラグメント対応部分を切り出す
                pdf.run_pdf([self._cmds['mat-select'],
                             '-t', global_dim,
                             '-l', global_dim,
                             '-b', global_dim + dim - 1,
                             '-r', global_dim + dim - 1,
                             dens_mat_path,
                             frg_dens_mat_path])

                # select された行列を対称行列に変換
                pdf.run_pdf([self._cmds['mat-symmetrize'],
                             frg_dens_mat_path,
                             frg_dens_mat_path])
                logger.debug("{header} density matrix for {fragment} was saved as {path}".format(
                    header=self.header,
                    fragment=frg_name,
                    path=frg_dens_mat_path))

                is_loadable = pdf.SymmetricMatrix.is_loadable(
                    frg_dens_mat_path)
                assert(is_loadable == True)
                (row, col) = pdf.SymmetricMatrix.get_size(frg_dens_mat_path)
                assert(row == dim)
                assert(row == col)

                # 対称行列パスをフラグメントに登録
                frg.set_density_matrix(frg_dens_mat_path)

                global_dim += dim

        logger.is_finished_pickup_density_matrix = True
        self.save()
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
        self.save()
        self.restore_cwd()

    # ------------------------------------------------------------------
    def pickup_QCLO_matrix(self, run_type='rks', force=False):
        if ((self.is_finished_pickup_LO == True) and
                (force == False)):
            logger.info('pickup LO has been finished.')
            return

        self.calc_lo(run_type)

        self.cd_work_dir('pickup lo')

        # post-SCF
        self._grouping_fragments()
        self._switch_fragments()

        # debug
        logger.debug("pickup_QCLO_matrix frame: ".format(self._name))
        pdfarc = self.get_pdfarchive()
        num_of_AOs = pdfarc.num_of_AOs
        num_of_MOs = pdfarc.num_of_MOs
        HOMO_level = pdfarc.get_HOMO_level('rks')  # option base 0
        logger.info('num of AOs: {}'.format(num_of_AOs))
        logger.info('num of MOs: {}'.format(num_of_MOs))
        logger.info('HOMO level: {}'.format(HOMO_level + 1))

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
        Clo_path = self.pdfparam.get_clo_mat_path()
        pdf.run_pdf(['component',
                     '-v',
                     '-S', 'CSC.mat',
                     '-c', Clo_path])

        # load CSC
        CSC = pdf.Matrix()
        CSC.load(CSC_path)

        logger.info("{header} make AO v.s. fragment table".format(
            header=self.header))
        AO_frg_tbl = self._get_AO_fragment_table(num_of_AOs)

        # pickup
        logger.info('{header} assign fragment: start: HOMO={homo}'.format(
            header=self.header, homo=HOMO_level))
        MO_fragment_assigned = {}
        for mo in range(HOMO_level + 1):
            frg_name = self._define_lo_fragment(
                mo, num_of_AOs, AO_frg_tbl, CSC)
            logger.info("{header} #{mo} MO -> fragment: '{frg_name}'".format(
                header=self.header, mo=mo, frg_name=frg_name))
            MO_fragment_assigned.setdefault(frg_name, [])
            MO_fragment_assigned[frg_name].append(mo)
        logger.info("{header} assign fragment: end".format(header=self.header))

        # assign report
        logger.info('==== assign report ====')
        for k, MOs in MO_fragment_assigned.items():
            logger.info("{header} fragment '{frag_name}' has {mo} MO(s)".format(
                header=self.header, frag_name=k, mo=len(MOs)))

        # フラグメントのC_LOを作成する
        logger.info("{header} create C_LO: start".format(header=self.header))
        Clo = pdf.Matrix()
        Clo.load(Clo_path)
        assert(num_of_AOs == Clo.rows)
        for frg_name, frg in self.fragments():
            frg_cols = len(MO_fragment_assigned.get(frg_name, []))
            logger.info("{header} fragment '{frg_name}': col={col}".format(
                header=self.header, frg_name=frg_name, col=frg_cols))
            if frg_cols == 0:
                logger.warning("{header} fragment '{frg_name}' has no colomns.".format(
                    header=self.header, frg_name=frg_name))
                # continue

            Clo_frg = pdf.Matrix(num_of_AOs, frg_cols)
            if frg_name in MO_fragment_assigned:
                for col, ref_col in enumerate(MO_fragment_assigned[frg_name]):
                    for row in range(num_of_AOs):
                        v = Clo.get(row, ref_col)
                        Clo_frg.set(row, col, v)

            Clo_path = 'Clo_{}.mat'.format(frg_name)
            logger.debug("{header} fragment C_LO save: {path}".format(
                header=self.header, path=Clo_path))
            Clo_frg.save(Clo_path)
            frg.set_LO_matrix(Clo_path, run_type)
        logger.info("{header} create C_LO: end".format(header=self.header))

        # trans C_LO to QCLO
        self._trans_LO2QCLO()

        # finish
        self.is_finished_pickup_LO = True
        self.save()
        self.restore_cwd()

    def _get_AO_fragment_table(self, num_of_AOs):
        '''
        AO v.s. fragment_name の辞書を返す
        '''
        frg_table = [None for x in range(num_of_AOs)]

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

        ranked_judge = sorted(judge.items(), key=lambda x: x[1], reverse=True)

        for rank, (k, v) in enumerate(ranked_judge):
            logger.info("{header} [{rank}] name:{name}, score:{score:.3f}".format(
                header=self.header,
                rank=rank + 1,
                name=k,
                score=v))

        high_score = ranked_judge[0][1]
        if high_score < 0.5:
            logger.warning("{header} 1st score is too small: {score}".format(
                header=self.header, score=high_score))

        return ranked_judge[0][0]

    def _trans_LO2QCLO(self):
        logger.info('trans LO at {}'.format(os.getcwd()))
        run_type = 'rks'
        F_path = self.pdfparam.get_f_mat_path(run_type)
        logger.info('F matrix: {}'.format(F_path))
        for frg_name, frg in self.fragments():
            C_QCLO_path = 'C_QCLO.{}.mat'.format(
                frg_name)  # output for each fragment
            frg_AO = frg.get_number_of_AOs()
            logger.info("{header} fragment '{name}' has {ao} AO(s)".format(
                header=self.header, name=frg_name, ao=frg_AO))
            if frg.get_number_of_AOs() != 0:
                Clo_path = frg.get_LO_matrix_path(run_type)
                assert(Clo_path != None)

                # calc (C_LO)dagger * F * C_LO => F'
                F_Clo_path = 'F_Clo.{}.mat'.format(frg_name)
                pdf.run_pdf([self._cmds['mat-mul'], '-v',
                             F_path,
                             Clo_path,
                             F_Clo_path])

                Clo_dagger_path = 'Clo_dagger.{}.mat'.format(frg_name)
                pdf.run_pdf([self._cmds['mat-transpose'], '-v',
                             Clo_path,
                             Clo_dagger_path])

                F_prime_path = 'Fprime.{}.mat'.format(frg_name)
                pdf.run_pdf([self._cmds['mat-mul'], '-v',
                             Clo_dagger_path,
                             F_Clo_path,
                             F_prime_path])

                pdf.run_pdf([self._cmds['mat-symmetrize'],
                             F_prime_path,
                             F_prime_path])

                # diagonal F'
                eigval_path = 'QCLO_eigval.{}.vtr'.format(frg_name)
                Cprime_path = 'Cprime.{}.mat'.format(frg_name)
                logger.info("diagonal F'")
                pdf.run_pdf([self._cmds['mat-diagonal'], '-v',
                             '-l', eigval_path,
                             '-x', Cprime_path,
                             F_prime_path])

                # AO基底に変換
                pdf.run_pdf([self._cmds['mat-mul'], '-v',
                             Clo_path,
                             Cprime_path,
                             C_QCLO_path])
            else:
                logger.info("{header} create empty QCLO matrix.".format(
                    header=self.header))
                empty_mat = pdf.Matrix()
                empty_mat.save(C_QCLO_path)

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
            logger.debug(
                "rearrangement of fragments is prohibited after calculation.")
            return

        if 'frame_molecule' in self._cache:
            self._cache.pop('frame_molecule')

        fragment_name = bridge.Utils.to_unicode(fragment_name)

        if isinstance(obj, QcFragment):
            fragment = QcFragment(obj)
            fragment.parent = self
            fragment.name = fragment_name
            logger.debug('[{my_name}] add fragment: name={fragment_name}'.format(
                my_name=self.name,
                fragment_name=fragment_name))
            self._fragments[fragment_name] = fragment
        elif isinstance(obj, QcFrame):
            logger.info(
                'begin to register frame molecule: for {}'.format(fragment_name))
            fragment = QcFragment()
            fragment.parent = self
            fragment.name = fragment_name
            for k, f in obj.fragments():
                if not f.margin:
                    logger.warn(
                        'add fragment: fragment={} for {}'.format(k, fragment_name))
                    fragment.set_group(k, f)
                else:
                    logger.warn(
                        'pass fragment: fragment={} is margin'.format(k))
            self._fragments[fragment_name] = fragment
            logger.info(
                'end of registration frame molecule: for {}'.format(fragment_name))
        else:
            raise

    # rearangement -----------------------------------------------------
    def _switch_fragments(self):
        '''
        fragmentsを入力用から出力用に切り替える

        処理内容:
        - 各fragmentの親を自分(self)にする
        '''
        logger.info("{header} switch fragment".format(header=self.header))
        output_fragments = OrderedDict()
        for frg_name, frg in self.fragments():
            logger.info("{header} fragment_name: {name}".format(
                header=self.header, name=frg_name))
            new_frg = QcFragment(frg, parent=self)
            assert(new_frg.parent.name == self.name)
            output_fragments[frg_name] = new_frg
        self._fragments = output_fragments

        # logger.info('merge subgroups')
        # for key, frg in self.fragments():
        #    frg.merge_subgroups()

        logger.info("{header} ---> switch".format(header=self.header))
        for frg_name, frg in self.fragments():
            logger.info("{header} {frg_name}: parent={parent_name}".format(
                header=self.header,
                frg_name=frg_name,
                parent_name=frg.parent.name))
        logger.info("{header} <---".format(header=self.header))

    def _grouping_fragments(self):
        logger.info("{header} grouping fragments".format(header=self.header))
        for frg_name, frg in self.fragments():
            frg.grouping_subfragments()

    # ==================================================================
    # coordinates
    # ==================================================================
    # outout XYZ -------------------------------------------------------
    def output_xyz(self, file_path):
        xyz = bridge.Xyz(self.frame_molecule)
        xyz.save(file_path)

    def check_bump_of_atoms(self):
        logger.info("{header} check bump of atoms".format(header=self.header))
        atom_list = self.frame_molecule.get_atom_list()
        num_of_atoms = len(atom_list)
        for i in range(num_of_atoms):
            xyz1 = atom_list[i].xyz
            for j in range(i):
                d = xyz1.distance_from(atom_list[j].xyz)
                if d < self.TOO_SMALL:
                    logger.warning("{header} atom[{i}][{atom_i}]({atom_i_path}) is near by atom[{j}][{atom_j}]({atom_j_path})".format(
                        header=self.header,
                        i=i, atom_i=str(atom_list[i]), atom_i_path=atom_list[i].path,
                        j=j, atom_j=str(atom_list[j]), atom_j_path=atom_list[j].path))
        logger.debug("{header} check_bump of atoms: end".format(
            header=self.header))

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

    # ==================================================================
    # debug
    # ==================================================================

    def _get_logger_header(self):
        """return header string for logger
        """
        header = "{name}>".format(name=self.name)
        return header
    header = property(_get_logger_header)
