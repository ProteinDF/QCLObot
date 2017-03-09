#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging
logger = logging.getLogger(__name__)

import pdfbridge
from .mdobject import MdObject
from .process import Process
from .utils import get_model, check_format_model, remove_WAT

class AmberObject(MdObject):
    """Operate Amber task

    : Example

    prepare molecule by AtomGroup(@pdfbridge)
    >>> pdb = pdfbridge.Pdb("./data/sample/1AKG.pdb")
    >>> models = pdb.get_atomgroup()
    >>> model = get_model(models)

    make working directory
    >>> md1 = AmberObject(name='md_test1')

    set model
    >>> md1.model = model

    optimization
    >>> md1.opt()
    
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    def _initialize(self):
        ''' initialize object
        '''
        super()._initialize() # called from the parents class

        self._AMBERHOME = os.environ.get('AMBERHOME', '')
        
        self._data['leap_sources'] = ['leaprc.protein.ff14SB',
                                      'leaprc.water.tip3p',
                                      'leaprc.gaff']
        self._data['leap_amberparams'] = []
        
    # ==================================================================
    # properties
    # ==================================================================
    def _get_model(self):
        return self._data.get('model', None)
    def _set_model(self, model):
        assert(isinstance(model, pdfbridge.AtomGroup))
        if check_format_model(model):
            self._data['model'] = pdfbridge.AtomGroup(model)
        else:
            logger.critical("not support the format; use model format")
            raise
    model = property(_get_model, _set_model)

    # solvation
    def _get_is_solvation(self):
        return self._data.get('is_solvation', False)
    def _set_is_solvation(self, yn):
        self._data['is_solvation'] = bool(yn)
    is_solvation = property(_get_is_solvation, _set_is_solvation)

    def _get_solvation_model(self):
        return self._data.get('solvation_model', None)
    def _set_solvation_model(self, model):
        self._data['solvation_model'] = str(model)
    solvation_model = property(_get_solvation_model, _set_solvation_model)

    # other amber parameters
    def set_param(self, param, value):
        self._data.setdefault("cntrl", {})
        param = str(param)
        value = str(value)
        self._data["cntrl"][param] = value
    
    # ==================================================================
    # properties: filepath
    # ==================================================================
    def _get_leapin_filepath(self):
        path = os.path.join(self.work_dir, 'leap.in')
        return path
    leapin_filepath = property(_get_leapin_filepath)

    def _get_leap_logfile_filepath(self):
        return os.path.join(self.work_dir,
                            self._data.get('logfile_filepath', 'leap.log'))
    leap_logfile_filepath = property(_get_leap_logfile_filepath)

    def _get_leap_sources(self):
        return self._data['leap_sources']
    leap_sources = property(_get_leap_sources)

    def _get_leap_amberparams(self):
        return self._data['leap_amberparams']
    leap_amberparams = property(_get_leap_amberparams)
    
    def _get_input_pdb_filepath(self):
        path = self._data.get('input_pdb_filepath', 'input.pdb')
        return path
    input_pdb_filepath = property(_get_input_pdb_filepath)

    def _get_output_pdb_filepath(self):
        path = self._data.get('output_pdb_filepath', 'output.pdb')
        return path
    output_pdb_filepath = property(_get_output_pdb_filepath)

    def _get_prmtop_filepath(self):
        path = self._data.get('prmtop_filepath', 'protein.prmtop')
        return path
    prmtop_filepath = property(_get_prmtop_filepath)

    def _get_inpcrd_filepath(self):
        path = self._data.get('inpcrd_filepath', 'protein.inpcrd')
        return path
    inpcrd_filepath = property(_get_inpcrd_filepath)
    
    def _get_initial_pdb_filepath(self):
        path = self._data.get('initial_pdb_filepath', 'initial.pdb')
        return path
    initial_pdb_filepath = property(_get_initial_pdb_filepath)

    def _get_mdin_filepath(self):
        path = self._data.get('mdin_filepath', 'md.in')
        return path
    mdin_filepath = property(_get_mdin_filepath)

    def _get_mdout_filepath(self):
        path = self._data.get('mdout_filepath', 'md.out')
        return path
    mdout_filepath = property(_get_mdout_filepath)

    def _get_restart_filepath(self):
        path = self._data.get('restart_filepath', "md.restart")
        return path
    restart_filepath = property(_get_restart_filepath)

    def _get_mdcrd_filepath(self):
        path = self._data.get("mdcrd_filepth", "mdcrd")
        return path
    mdcrd_filepath = property(_get_mdcrd_filepath)

    def _get_mdvel_filepath(self):
        path = self._data.get("mdvel_filepath", "mdvel")
        return path
    mdvel_filepath = property(_get_mdvel_filepath)

    def _get_mden_filepath(self):
        path = self._data.get("mden_filepath", "mden")
        return path
    mden_filepath = property(_get_mden_filepath)

    def _get_final_pdb_filepath(self):
        path = self._data.get('final_pdb_filepath', 'final.pdb')
        return path
    final_pdb_filepath = property(_get_final_pdb_filepath)
    
    # ==================================================================
    # internal properties
    # ==================================================================
    def _set_chainres_match_table(self,
                                  amb_chain, amb_res,
                                  orig_chain, orig_res):
        self._data.setdefault("cr_tbl", {})
        key = (amb_chain, amb_res)
        self._data["cr_tbl"][key] = (orig_chain, orig_res)

    def _get_chainres_match_table(self, amb_chain, amb_res):
        self._data.setdefault("cr_tbl", {})
        key = (amb_chain, amb_res)
        return self._data["cr_tbl"].get(key, (None, None))
    
    # ==================================================================
    # optimization
    # ==================================================================
    def opt(self):
        model = pdfbridge.AtomGroup(self.model)
        model = self._remove_water(model)
        self._prepare_input_pdb(model)
        self._prepare_leapin()
        self._do_leap()

        self.set_param("ntx", 1)

        self.set_param("ntpr", 200)

        self.set_param("ntc", 1)

        self.set_param("ntf", 1)
        self.set_param("ntb", 0)
        self.set_param("cut", 9999.999)
        self.set_param("igb", 1)

        self._prepare_mdin()

        self._do_sander()
        self._restrt2pdb()
        
    def _remove_water(self, model):
        return remove_WAT(model)

        
    def _prepare_input_pdb(self, model):
        self.cd_workdir("prepare input pdb file")

        #mode = ''
        mode = 'amber'

        models = pdfbridge.AtomGroup()
        models.set_group('model_1', model)
        pdb = pdfbridge.Pdb(mode=mode)
        pdb.set_by_atomgroup(models)

        with open(self.input_pdb_filepath, "w") as f:
            f.write(str(pdb))

        self.restore_cwd()

    # ==================================================================
    # dynamics
    # ==================================================================
    def md(self,
           steps = 1, dt = 0.002):
        """ compute molecular dynamics

        steps: Number of MD-steps to be performed.
        dt:    The time step (psec).
        """
        model = pdfbridge.AtomGroup(self.model)
        model = self._remove_water(model)
        self._prepare_input_pdb(model)
        self.solvation_model = "cap"
        
        self._prepare_leapin()
        self._do_leap()

        #self.set_param("imin", 0)
        self.set_param("irest", 0)

        self.set_param("ntpr", 5000)
        self.set_param("ntwx", 500)
        self.set_param("ntwe", 500)

        self.set_param("nstlim", steps) # 500000
        self.set_param("dt", dt)

        self.set_param("ntt", 1)
        self.set_param("temp0", 300.0)
        self.set_param("tempi", 0.0)
        self.set_param("tautp", 0.1)

        # SHAKE bond length constraints
        self.set_param("ntc", 2)

        self.set_param("ntf", 2)
        self.set_param("ntb", 0)
        self.set_param("cut", 10.0)

        self._prepare_mdin(opt=False)

        self._do_sander()
        self._restrt2pdb()
        
        
    # ==================================================================
    # leap setting
    # ==================================================================
    def _prepare_leapin(self):
        """make input file for xleap command
        """
        self.cd_workdir("prepare leapin")
        
        leapin_contents = ''
        leapin_contents += 'logFile {logfile_filepath}\n'.format(logfile_filepath=self.leap_logfile_filepath)
        leapin_contents += self.__get_leap_source_lines()
        leapin_contents += self.__get_leap_amberparams_lines()
        leapin_contents += 'protein = loadPdb {pdb_file}\n'.format(pdb_file=self.input_pdb_filepath)
        leapin_contents += 'proteinBox = copy protein\n'
        leapin_contents += self.__get_leap_ssbond_lines()
        leapin_contents += self.__get_solvation(solute='proteinBox', solvent='TIP3PBOX')
        leapin_contents += 'saveAmberParm proteinBox {prmtop} {inpcrd}\n'.format(prmtop=self.prmtop_filepath,
                                                                                 inpcrd=self.inpcrd_filepath)
        leapin_contents += 'savePdb proteinBox {pdb_file}\n'.format(pdb_file=self.initial_pdb_filepath)
        leapin_contents += 'quit\n'
        
        logger.debug('save leap inputfile: {}'.format(self.leapin_filepath))
        with open(self.leapin_filepath, 'w') as f:
            f.write(leapin_contents)

        self.restore_cwd()
            
    def __get_leap_source_lines(self):
        """return source command lines for leap.in
        """
        answer = ''
        for ff in self.leap_sources:
            answer += "source {ff}\n".format(ff=ff)
        return answer

    def __get_leap_amberparams_lines(self):
        """return loadamberparams and loadamberprep lines for leap.in
        """
        answer = ''
        for lig in self.leap_amberparams:
            answer += 'loadAmberParams {lig}.frcmod\n'.format(lig=lig)
            answer += 'loadAmberPrep {lig}.prep\n'.format(lig=lig)
        return answer

    def __get_leap_ssbond_lines(self):
        """return bond lines
        """
        answer = ''

        return answer

    def __get_solvation(self, solute, solvent):
        """return solvation leap code.
        """
        answer = ''

        if self.solvation_model != None:
            if self.solvation_model.upper() == 'CAP':
                center = self.model.center()
                cap_center = "{{ {x:.3f} {y:.3f} {z:.3f} }}".format(x=center.x,
                                                                    y=center.y,
                                                                    z=center.z)
                (box_min, box_max) = self.model.box()
                distance = 0.0
                distance = max(distance, center.distance_from(pdfbridge.Position(box_min.x, box_min.y, box_min.z)))
                distance = max(distance, center.distance_from(pdfbridge.Position(box_max.x, box_min.y, box_min.z)))
                distance = max(distance, center.distance_from(pdfbridge.Position(box_max.x, box_max.y, box_min.z)))
                distance = max(distance, center.distance_from(pdfbridge.Position(box_max.x, box_min.y, box_max.z)))
                distance = max(distance, center.distance_from(pdfbridge.Position(box_max.x, box_max.y, box_max.z)))
                distance = max(distance, center.distance_from(pdfbridge.Position(box_min.x, box_max.y, box_min.z)))
                distance = max(distance, center.distance_from(pdfbridge.Position(box_min.x, box_max.y, box_max.z)))
                distance = max(distance, center.distance_from(pdfbridge.Position(box_min.x, box_min.y, box_max.z)))
                cap_radius = max(distance + 10.0, 30.0)
                answer += "solvateCap {solute} {solvent} {position} {radius} {closeness}\n".format(
                    solute=solute,
                    solvent=solvent,
                    position=cap_center,
                    radius=cap_radius,
                    closeness="")
            else:
                logger.warning('solvation model is not understood.: {}'.format(self.solvation_model))
        
        return answer
    
    # ==================================================================
    # do leap
    # ==================================================================
    def _do_leap(self):
        self.cd_workdir("do leap")
        
        p = Process()
        leap_cmd = os.path.join(self._AMBERHOME, 'bin', 'tleap')
        cmd = "{leap_cmd} -s -f {leapin}".format(leap_cmd=leap_cmd,
                                                 leapin=self.leapin_filepath)
        p.cmd(cmd)
        return_code = p.commit(stdout_through=False,
                               stderr_through=False)

        self.restore_cwd()
        return return_code

    # ==================================================================
    # prepare mdin
    # ==================================================================
    def _prepare_mdin(self,
                      opt=True):
        self.cd_workdir("prepare md inputfile")

        imin = 1 if opt else 0
        basic_input_options = "imin={imin}".format(
            imin=imin
        )

        ntx = 1
        format_of_input = "ntx={ntx}".format(
            ntx=ntx
        )

        minimization_contents = ""
        lmod_contents = ""
        if opt:
            # maxcyc: The maximum number of cycles of minimization. Default = 1.
            # ncyc: the method of minimization will be switched from steepest
            #   descent to conjugate gradient after NCYC cycles. (default:10)
            # ntmin: Flag for the minimization (Amber)
            #   = 0 Full conjugate gradient minimization. The first 4 cycles are steepest descent at
            #       the start of the run and after every nonbonded pairlist update.
            #   = 1 For NCYC cycles the steepest descent method is used then conjugate gradient
            #       is switched on (default).
            #   = 2 Only the steepest descent method is used.
            #   = 3 The XMIN method is used.
            #   = 4 The LMOD method is used.
            # dx0: The initial step length. (default 0.01)
            # drms: The convergence criterion for the energy gradient. (default 1.0E-4 kcal/mol A)
            # xmin_method:
            #   "PRCG" = Polak-Ribiere Conjugate Gradient,
            #   "LBFGS" = Limited-memory Broyden-Fletcher-Goldfarb-Shanno (default)
            #   "TNCG" = Optionally LBFGS-preconditioned Truncated Newton Conjugate Gradient.
            maxcyc = 50000
            ncyc = 0
            ntmin = 3
            drms = 0.001
            # dx を指定するとNG
            minimization_contents = "imin=1, maxcyc={maxcyc}, ncyc={ncyc}, ntmin={ntmin}, drms={drms},".format(
                maxcyc=maxcyc,
                ncyc=ncyc,
                ntmin=ntmin,
                drms=drms
            )
            lmod_contents = """
        &lmod
          xmin_method="LBFGS"
        &end
            """
        else:
            # molecular dynamics
            # nstlim: Number of MD-steps to be performed. Default 1.
            minimization_contents = "imin=0, "

            
        # ntx: Option to read the coordinates from the “inpcrd” file.
        #  = 1 X is read formatted with no initial velocity information. Default.
        #  = 2 X is read unformatted with no initial velocity information.
        #
        # ntpr: Print the progress of the minimization every ntpr steps; default is 10.
        # 
        # **** Potential function parameters ****
        # ntb: periodic boundaries or not
        #  = 0 noperiodicityisappliedandPMEisoff(defaultwhenigb>0)
        #  = 1 constantvolume(defaultwhenigbandntpareboth0,whicharetheirdefaults)
        #  = 2 constant pressure (default when ntp > 0)
        mdin_contents = """
        # 
        &cntrl
          {minimization_contents}
          {cntrl_contents}
        &end
        """

        cntrl_contents = ""
        for param, value in self._data["cntrl"].items():
            cntrl_contents += "{param}={value}, ".format(
                param=param,
                value=value)
        
        mdin_contents = mdin_contents.format(
            minimization_contents=minimization_contents,
            cntrl_contents=cntrl_contents)
        mdin_contents += lmod_contents
        mdin_contents = mdin_contents.lstrip('\n')
        mdin_contents = pdfbridge.Utils.unindent_block(mdin_contents)

        with open(self.mdin_filepath, 'w') as f:
            f.write(mdin_contents)

        self.restore_cwd()

#        mdin_contents = """
#        &cntrl
#          ibelly=1,
#          bellymask='${belly_mask}'
#        &end        
#        """
        
        
    # ==================================================================
    # do sander
    # ==================================================================
    def _do_sander(self):
        self.cd_workdir("do sander")

        p = Process()
        sander_cmd = os.path.join(self._AMBERHOME, "bin", "sander")
        cmd = "{sander_cmd} -i {mdin} -o {mdout} -c {inpcrd} -p {prmtop} -r {restart} -x {mdcrd} -v {mdvel} -e {mden} ".format(
            sander_cmd=sander_cmd,
            mdin=self.mdin_filepath,
            mdout=self.mdout_filepath,
            inpcrd=self.inpcrd_filepath,
            prmtop=self.prmtop_filepath,
            restart=self.restart_filepath,
            mdcrd=self.mdcrd_filepath,
            mdvel=self.mdvel_filepath,
            mden=self.mden_filepath)

        logger.debug("run command: {}".format(cmd))
        p.cmd(cmd)
        return_code = p.commit()
        
        self.restore_cwd()
        return return_code


    # ==================================================================
    # get pdb
    # ==================================================================
    def _restrt2pdb(self):
        self._make_matching_table()
        
        self.cd_workdir("restrt2pdb")

        p = Process()
        ambpdb_cmd = os.path.join(self._AMBERHOME, "bin", "ambpdb")
        cmd = "{ambpdb_cmd} -p {prmtop} -c {restrt}".format(ambpdb_cmd=ambpdb_cmd,
                                                            prmtop=self.prmtop_filepath,
                                                            restrt=self.restart_filepath)
        
        logger.debug("run command: {}".format(cmd))
        p.cmd(cmd)
        return_code = p.commit(self.final_pdb_filepath,
                               stdout_through=False,
                               stderr_through=False)

        amb_pdb = pdfbridge.Pdb(self.final_pdb_filepath)
        amb_models = amb_pdb.get_atomgroup()
        amb_model = get_model(amb_models)

        orig_model = self._rename_amb2orig(amb_model)
        orig_models = pdfbridge.AtomGroup()
        orig_models.set_group(1, orig_model)
        orig_pdb = pdfbridge.Pdb()
        orig_pdb.set_by_atomgroup(orig_models)
        with open(self.output_pdb_filepath, "w") as f:
            f.write(str(orig_pdb))
        
        self.restore_cwd()
        return return_code
    
    def _rename_amb2orig(self, amb_model):
        # rename model, chain, res_id
        answer = pdfbridge.AtomGroup()
        for amb_chain_id, amb_chain in amb_model.groups():
            for amb_res_id, amb_res in amb_chain.groups():
                orig_chain_id, orig_res_id = self._get_chainres_match_table(amb_chain_id, amb_res_id)
                if (orig_chain_id == None) or (orig_res_id == None):
                    orig_chain_id = amb_chain_id
                    orig_res_id = amb_res_id
                    
                if answer.has_groupkey(orig_chain_id) != True:
                        orig_chain = pdfbridge.AtomGroup(name=amb_chain.name)
                        answer.set_group(orig_chain_id, orig_chain)
                answer[orig_chain_id].set_group(orig_res_id, amb_res)
        
        return answer
    
    # ==================================================================
    # make matching table
    # ==================================================================
    def _make_matching_table(self):
        """
        make matching table between original and leap model.

        入力したモデルとleapで処理したモデルで残基番号などが変化するので、
        対応表を作成する
        """
        self.cd_workdir("match models")

        # make model_amber
        amb_pdb = pdfbridge.Pdb(self.initial_pdb_filepath)
        amb_models = amb_pdb.get_atomgroup()
        amb_model = get_model(amb_models)
        assert(check_format_model(amb_model))

        # check wether 'leap' adds any atoms or not.
        if self.model.get_number_of_all_atoms() == amb_model.get_number_of_all_atoms():
            logger.debug('Check number of atoms. It looks good.')
        elif self.model.get_number_of_all_atoms() < amb_model.get_number_of_all_atoms():
            logger.info('some atoms are added via the leap treatment.: original={} new={}'.format(
                self.model.get_number_of_all_atoms(), amb_model.get_number_of_all_atoms()))
        else:
            logger.infog('some atoms are deleteed via the leap treatment.: original={} new={}'.format(
                self.model.get_number_of_all_atoms(), amb_model.get_number_of_all_atoms()))
            
        # 入力モデルとleap後のモデルの座標は変わらないので、
        # 座標をもとに対応表を作成する

        # 座標が等しい原子をmodelから探す
        def find_atom_path(model, pos):
            assert(isinstance(model, pdfbridge.AtomGroup))
            assert(isinstance(pos, pdfbridge.Position))

            NEAR_DISTANCE = 0.1
            range_selecter = pdfbridge.Select_Range(pos, NEAR_DISTANCE)
            selection = model.select(range_selecter)
            path_list = selection.get_path_list()
            
            answer = None
            if len(path_list) > 0:
                answer = path_list.pop(0)
                
            return answer
                        
        # make original table
        # self._input_amber_model_match_table = []
        for chain_id, chain in amb_model.groups():
            for res_id, res in chain.groups():
                for atom_id, atom in res.atoms():
                    amb_path = atom.path
                    orig_path = find_atom_path(self.model, atom.xyz)

                    if orig_path == None:
                        pass
                    else:
                        (orig_model_id, orig_chain_id, orig_resid, orig_atom_id) = pdfbridge.AtomGroup.divide_path(orig_path)
                        (amb_model_id, amb_chain_id, amb_resid, amb_atom_id) = pdfbridge.AtomGroup.divide_path(amb_path)
                        self._set_chainres_match_table(amb_chain_id, amb_resid, orig_chain_id, orig_resid)

        self._show_match_table()
        
        self.restore_cwd()
        return 0

    
    def _show_match_table(self):
        ''' output matching table
        '''
        logger.debug(">>>> matching table")
        for (amb_chain_id, amb_resid), (orig_chain_id, orig_resid) in self._data["cr_tbl"].items():
            logger.debug("amb: {}/{} <-> orig: {}/{}".format(amb_chain_id, amb_resid, orig_chain_id, orig_resid))

            
if __name__ == '__main__':
    logger.setLevel(logging.DEBUG)
    import doctest
    doctest.testmod()


