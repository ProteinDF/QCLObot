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
import shelve


try:
    import msgpack
except:
    import msgpack_pure as msgpack

import pdfbridge as bridge
import pdfpytools as pdf
import qclobot as qclo

class QcFrame(object):
    _pdfparam_filename = 'pdfparam.mpac'
    _db_filename = 'pdfresults.db'

    # ------------------------------------------------------------------
    def __init__(self, name, *args, **kwargs):
        '''
        create QcFrame object.
        
        name: name of the frame molecule
        '''
        self._logger = logging.getLogger(__name__)

        # mandatory parameter
        self._name = name
        self._input_fragments = OrderedDict()
        self._output_fragments = OrderedDict()
        self._state = {} # 状態保持
        self._select_fragments()
        
        self._prepare_work_dir()
        self._load()

        # cache data
        self._frame_molecule = None
        self._pdfparam = None

        if ((len(args) > 0) and isinstance(args[0], QcFrame)):
            self._copy_constructer(args[0])
            
    # copy constructer
    def _copy_constructer(self, rhs):
        self._name = rhs._name
        self._pdfparam = rhs._pdfparam

    def __del__(self):
        self._save()
        pass
        
    def _load(self):
        return
        data = shelve.open(os.path.join(self.work_dir, 'qcframe'))
        self._fragments = data.get('fragments', {})
        data.close()
        
    def _save(self):
        return 
        data = shelve.open(os.path.join(self.work_dir, 'qcframe'))
        for k, v in self._fragments.items():
            print(k, v)
        data['fragments'] = self._fragments
        data.close()
        
    # pdfparam ---------------------------------------------------------
    def _get_pdfparam(self):
        '''
        pdfparamオブジェクトを返す
        '''
        pdfparam_path = os.path.abspath(os.path.join(self.work_dir,
                                                     self._pdfparam_filename))

        if self._pdfparam is None:
            pdfparam = None
            if os.path.exists(pdfparam_path):
                mpac_file = open(pdfparam_path, 'rb')
                mpac_data = msgpack.unpackb(mpac_file.read())
                mpac_data = bridge.Utils.byte2str(mpac_data)
                mpac_file.close()
                self._pdfparam = pdf.PdfParam(mpac_data)
            else:
                pdfsim = pdf.PdfSim()
                self._pdfparam = pdf.get_default_pdfparam()
                #self._pdfparam.save(pdfparam_path)
            
        return self._pdfparam
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
        if self._frame_molecule == None:
            self._logger.info('create frame molecule coordinates.')
            frame_molecule = bridge.AtomGroup()
            for frg_name, frg in self._input_fragments.items():
                self._logger.info('fragment name={}: {} atoms'.format(frg_name,
                                                                      frg.get_number_of_all_atoms()))
                frame_molecule[frg_name] = frg.get_AtomGroup()
            self._frame_molecule = frame_molecule
            self._logger.info('')
        return self._frame_molecule

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
            self._logger.info('make workdir: {}'.format(self.work_dir))
            os.mkdir(self.work_dir)
        else:
            self._logger.debug('already exist: {}'.format(self.work_dir))

    def cd_work_dir(self, job_name=''):
        '''
        作業ディレクトリをオブジェクトのwork_dirに移動する
        '''
        self._logger.info('=' * 20)
        self._logger.info('>>>> {job_name}@{frame_name}'.format(job_name = job_name,
                                                                frame_name = self.name))
        self._logger.info('=' * 20)
        self._prev_dir = os.path.abspath(os.curdir)
        os.chdir(self.work_dir)

    def restore_cwd(self):
        '''
        self.cd_work_dir() 以前のディレクトリに戻す
        '''
        os.chdir(self._prev_dir)
        self._logger.info('<<<<\n')

    def _check_path(self, path):
        if not os.path.exists(path):
            self._logger.warn('NOT FOUND: {}'.format(path))
        
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
        
    # ==================================================================
    # GUESS
    # ==================================================================
    # guess density ----------------------------------------------------
    def guess_density(self, run_type ='rks'):
        self.cd_work_dir('guess_density')
        guess_density_matrix_path = 'guess.density.{}.mat'.format(run_type)

        # 既存のデータを消去する
        if os.path.exists(guess_density_matrix_path):
            os.remove(guess_density_matrix_path)
        
        pdfsim = pdf.PdfSim()
        pdfsim.setup()
        
        for frg_name, frg in self.fragments():
            self._logger.info('fragment name={}: {} atoms'.format(frg_name,
                                                                  frg.get_number_of_all_atoms()))
            if frg.qc_parent == None:
                self._logger.warn('guess_density(): qc_parent == None. frg_name={}'.format(frg_name))

            frg_guess_density_matrix_path = frg.get_guess_density_matrix(run_type)

            self._logger.debug('guess_density() [{}@{}] ext: {} from {}'.format(
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
                self._logger.warn('not found: frg.guess.dens.mat={}'.format(frg_guess_density_matrix_path))

        self.pdfparam.guess = 'density_matrix'
        self._logger.info('initial guess (density matrix) created at {}'.format(guess_density_matrix_path))

        # check
        self._check_path(guess_density_matrix_path)
        
        self.restore_cwd()

        
    def guess_QCLO(self, run_type='rks', isCalcOrthogonalize = False):
        self.cd_work_dir('guess_QCLO')

        guess_QCLO_matrix_path = 'guess.QCLO.{}.mat'.format(run_type)
        if os.path.exists(guess_QCLO_matrix_path):
            os.remove(guess_QCLO_matrix_path)
        
        num_of_AOs = 0
        for frg_name, frg in self.fragments():
            self._logger.debug('guess QCLO: frg_name={}, parent={}'.format(frg_name, frg.qc_parent.name))
            frg_QCLO_matrix_path = frg.prepare_QCLO_matrix(run_type, self)
            if os.path.exists(frg_QCLO_matrix_path):
                pdf.run_pdf(['mat-ext', '-c',
                             guess_QCLO_matrix_path,
                             frg_QCLO_matrix_path,
                             guess_QCLO_matrix_path])
            else:
                self._logger.warn('The QCLO of the subgroup, {}, was not created.'.format(frg_name))

        # orthogonalize
        guess_path = 'guess.lcao.{}.mat'.format(run_type)
        if isCalcOrthogonalize:
            self._logger.info('orthogonalize')
            Xinv_path = self.pdfparam.get_Xinv_mat_path()

            self._check_path(guess_QCLO_matrix_path)
            pdf.run_pdf(['mat-mul', '-v',
                         Xinv_path,
                         guess_QCLO_matrix_path,
                         guess_path])
        else:
            shutil.copy(guess_QCLO_matrix_path, guess_path)

        self.pdfparam.guess = 'lcao'
        self._logger.info('guess LCAO matrix created: {}'.format(guess_path))

        # check
        self._check_path(guess_QCLO_matrix_path)

        self.create_occupation_file(run_type)
        
        self.restore_cwd()

    def create_occupation_file(self, run_type='rks'):
        self.cd_work_dir('create occ')

        occ_level = -1
        electrons_per_orb = 0.0
        run_type = run_type.upper()
        if run_type == 'RKS':
            occ_level = int((self.pdfparam.num_of_electrons / 2.0))
            electrons_per_orb = 2.0
        else:
            self._logger.critical('NOT supported. run_type={}'.format(run_type))

        num_of_MOs = self.pdfparam.num_of_MOs
        occ_vtr = pdf.Vector(num_of_MOs)
        for i in range(occ_level):
            occ_vtr.set(i, electrons_per_orb)

        occ_vtr_path = 'guess.occ.{}.vtr'.format(run_type.lower())
        occ_vtr.save(occ_vtr_path)
        self._check_path(occ_vtr_path)
            
        self.restore_cwd()
        
    # ==================================================================
    # CALC
    # ==================================================================
    # ------------------------------------------------------------------
    def calc_preSCF(self, dry_run=False):
        '''
        '''
        self.cd_work_dir('calc SP')

        pdfsim = pdf.PdfSim()

        for frg_name, frg in self.fragments():
            frg.set_basisset(self.pdfparam)
        self.pdfparam.molecule = self.frame_molecule

        self.pdfparam.step_control = 'integral'
        pdfsim.sp(self.pdfparam,
                  workdir = self.work_dir,
                  db_path = self.db_path,
                  dry_run = dry_run)

        self._pdfparam = None # cache clear
        self.is_finished_prescf = True
        
        self.restore_cwd()
    
    # sp ---------------------------------------------------------------
    def calc_sp(self, dry_run=False):
        '''
        calculate single point energy
        '''
        self.cd_work_dir('calc SP')

        pdfsim = pdf.PdfSim()

        if self.is_finished_scf != True:
            for frg_name, frg in self.fragments():
                frg.set_basisset(self.pdfparam)
            self.pdfparam.molecule = self.frame_molecule

        #self.output_xyz("{}/model.xyz".format(self.name))

        step_control = 'guess scf'
        if self.is_finished_prescf != True:
            step_control = 'integral ' + step_control
        self.pdfparam.step_control = step_control
        
        pdfsim.sp(self.pdfparam,
                  workdir = self.work_dir,
                  db_path = self.db_path,
                  dry_run = dry_run)

        self._pdfparam = None # cache clear
        self.is_finished_prescf = True
        self.is_finished_scf = True
        self._switch_fragments()
        
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
        self.cd_work_dir('pickup density matrix')
        
        dens_mat_path = self.pdfparam.get_density_matrix_path(runtype=runtype)
        self._logger.info('dens.mat={}'.format(dens_mat_path))
        
        global_dim = 0
        for frg_name, frg in self.fragments():
            dim = frg.get_number_of_AOs()
            if dim > 0:
                frg_dens_mat_path = 'Ppq.{}.{}.mat'.format(runtype, frg_name)
                self._logger.info(
                    'density matrix select frg:{} in frame:{} was saved as {} ({} -> {})'.format(
                        frg_name,
                        self.name,
                        frg_dens_mat_path,
                        global_dim,
                        global_dim +dim -1))

                # フラグメント対応部分を切り出す
                pdf.run_pdf(['mat-select', '-v',
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

                self._logger.info(
                    'density matrix of the frg:{} in frame:{} was saved as {}'.format(
                        frg_name,
                        self.name,
                        frg_dens_mat_path))
                global_dim += dim

        self.restore_cwd()
            
    # ------------------------------------------------------------------
    def calc_lo(self):
        self.cd_work_dir('calc lo')

        self._logger.info('start lo calculation.')
        pdf.run_pdf('lo')

        self.restore_cwd()

    # ------------------------------------------------------------------
    def pickup_lo(self):
        self.cd_work_dir('pickup lo')
        
        # debug
        pdfarc = self.get_pdfarchive()
        num_of_AOs = pdfarc.num_of_AOs
        num_of_MOs = pdfarc.num_of_MOs
        HOMO_level = pdfarc.get_HOMO_level('rks') # option base 0
        self._logger.info('num of AOs: {}'.format(num_of_AOs))
        self._logger.info('num of MOs: {}'.format(num_of_MOs))
        self._logger.info('HOMO level: {}'.format(HOMO_level +1)) 

        self._logger.info('fragment information:')
        for frg_name, frg in self.fragments():
            frg_AOs = frg.get_number_of_AOs()
            self._logger.info('fragment name:[{}] AOs={}'.format(frg_name, frg_AOs))
        self._logger.info('')
            
        # calc S*C
        lo_satisfied = self.pdfparam.lo_satisfied
        if lo_satisfied != True:
            self._logger.warn('lo_satisfied: {}'.format(lo_satisfied))
        lo_iterations = self.pdfparam.lo_num_of_iterations
        self._logger.info('lo iterations: {}'.format(lo_iterations))

        self._logger.info('calc S*C')
        CSC_path = 'CSC.mat'
        Clo_path = self.pdfparam.get_clomat_path()
        pdf.run_pdf(['component',
                     '-v',
                     '-S', 'CSC.mat',
                     '-c', Clo_path])

        # load CSC
        CSC = pdf.Matrix()
        CSC.load(CSC_path)

        self._logger.info('make AO v.s. fragment table')
        AO_frg_tbl = self._get_AO_fragment_table(num_of_AOs)
        
        # pickup
        self._logger.info('assign fragment: start: HOMO={}'.format(HOMO_level))
        MO_fragment_assigned = {}        
        for mo in range(HOMO_level +1):
            frg_name = self._define_lo_fragment(mo, num_of_AOs, AO_frg_tbl, CSC)
            self._logger.info('{}-th MO -> fragment: [{}]'.format(mo, frg_name))
            MO_fragment_assigned.setdefault(frg_name, [])
            MO_fragment_assigned[frg_name].append(mo)
        self._logger.info('assign fragment: end')

        # assign report
        self._logger.info('==== assign report ====')
        for k, MOs in MO_fragment_assigned.items():
            self._logger.info('fragment[{}] # of MOs={}'.format(k, len(MOs)))
        self._logger.info('')
            
        # フラグメントのC_LOを作成する
        self._logger.info('create C_LO: start')
        Clo = pdf.Matrix()
        Clo.load(Clo_path)
        assert(num_of_AOs == Clo.rows)
        for frg_name, frg in self.fragments():
            frg_cols = len(MO_fragment_assigned.get(frg_name, []))
            if frg_cols == 0:
                continue

            Clo_frg = pdf.Matrix(num_of_AOs, frg_cols)
            for col, ref_col in enumerate(MO_fragment_assigned[frg_name]):
                for row in range(num_of_AOs):
                    v = Clo.get(row, ref_col)
                    Clo_frg.set(row, col, v)

            Clo_path = 'Clo_{}.mat'.format(frg_name)
            self._logger.info('fragment C_LO save: {}'.format(Clo_path))
            Clo_frg.save(Clo_path)
            frg.Clo_path = Clo_path
        self._logger.info('create C_LO: end')
            
        # trans C_LO
        self._trans_LO()

        # finish
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
            self._logger.info('[{rank}] name:{name}, score:{score:.3f}'.format(
                rank=rank +1,
                name=k,
                score=v))

        high_score = ranked_judge[0][1]
        if high_score < 0.5:
            self._logger.warn('1st score is too small: {}'.format(high_score))
            
        return ranked_judge[0][0]

    def _trans_LO(self):
        self._logger.info('trans LO at {}'.format(os.getcwd()))
        run_type = 'rks'
        F_path = self.pdfparam.get_Fmat_path(run_type)
        self._logger.info('F matrix: {}'.format(F_path))
        for frg_name, frg in self.fragments():
            Clo_path = frg.Clo_path
            if Clo_path is None:
                continue

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
            self._logger.info("diagonal F'")
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

            self._logger.info('C_QCLO saved: {}'.format(C_QCLO_path))
        
        
    #  =================================================================
    #  for fragments
    #  =================================================================
    def fragments(self):
        '''
        フラグメントの名前とオブジェクトを返すイテレータ
        '''
        for k in self._fragments.keys():
            yield(k, self._fragments[k])
    
    # operator[] -------------------------------------------------------
    def __getitem__(self, fragment_name):
        '''
        出力用[]演算子
        '''
        fragment_name = str(fragment_name)
        return self._fragments.get(fragment_name, None)

    def __setitem__(self, fragment_name, fragment):
        '''
        入力用[]演算子

        計算前であれば代入可能(つまりモデリング中)であるが、
        計算後は代入できない
        '''
        if self.is_finished_scf:
            self._logger.warn('operator[] called after simulation.')
            return

        self._frame_molecule = None # clear molecule cache

        fragment_name = str(fragment_name)
        fragment = qclo.QcFragment(fragment)
        fragment.name = fragment_name
        if fragment.qc_parent == None:
            fragment.qc_parent = self
        self._input_fragments[fragment_name] = fragment

    # select fragments -------------------------------------------------
    def _select_fragments(self):
        '''
        出力用フラグメントを状況に応じて切り替える
        '''
        if self.is_finished_scf:
            self._fragments = self._output_fragments
        else:
            self._fragments = self._input_fragments
        
    # rearangement -----------------------------------------------------
    def _switch_fragments(self):
        '''
        fragmentsを入力用から出力用に切り替える

        処理内容:
        - 各fragmentの親を自分(self)にする
        '''
        self._logger.info('switch fragment.')
        output_fragments = OrderedDict()
        for frg_name, frg in self.fragments():
            new_frg = qclo.QcFragment(frg, qc_parent=self)
            assert(new_frg.qc_parent.name == self.name)
            output_fragments[frg_name] = new_frg
        self._input_fragments = self._fragments
        self._output_fragments = output_fragments
        self._fragments = self._output_fragments

        self._logger.info('---> switch ')
        for frg_name, frg in self.fragments():
            self._logger.info('{}: parent={}'.format(frg_name,
                                                     frg.qc_parent.name))
        self._logger.info('<---')
        

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

