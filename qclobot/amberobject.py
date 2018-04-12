#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging
logger = logging.getLogger(__name__)
import time

import pdfbridge
from .mdobject import MdObject
from .process import Process
from .utils import get_model, check_format_model, find_max_chain_id, remove_WAT

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
        super(AmberObject, self).__init__(*args, **kwargs)

    def _initialize(self):
        ''' initialize object
        '''
        super(AmberObject, self)._initialize() # called from the parents class

        self._AMBERHOME = os.environ.get('AMBERHOME', '')
        if len(self._AMBERHOME) == 0:
            logger.warn("The environment variable `AMBERHOME` is empty.")

        self._data['leap_sources'] = ['leaprc.protein.ff14SB',
                                      'leaprc.water.tip3p',
                                      'leaprc.gaff']
        self._data['leap_amberparams'] = []

    # ==================================================================
    # properties
    # ==================================================================
    # solvation
    #def _get_is_solvation(self):
    #    return self._data.get('is_solvation', False)
    #def _set_is_solvation(self, yn):
    #    self._data['is_solvation'] = bool(yn)
    #is_solvation = property(_get_is_solvation, _set_is_solvation)

    def _get_solvation_method(self):
        return self._data.get("solvation_method", None)
    def _set_solvation_method(self, method):
        self._data["solvation_method"] = str(method).lower()
    solvation_method = property(_get_solvation_method, _set_solvation_method)

    def _get_solvation_model(self):
        return self._data.get('solvation_model', "TIP3PBOX")
    def _set_solvation_model(self, model):
        self._data['solvation_model'] = str(model)
    solvation_model = property(_get_solvation_model, _set_solvation_model)


    # belly type dynamics
    def _get_use_belly(self):
        return self._data.get("use_belly", False)
    def _set_use_belly(self, yn):
        self._data["use_belly"] = bool(yn)
    use_belly = property(_get_use_belly, _set_use_belly)

    def _get_bellymask_WAT(self):
        return self._data.get("bellymask_WAT", False)
    def _set_bellymask_WAT(self, yn):
        self._data["bellymask_WAT"] = bool(yn)
    bellymask_WAT = property(_get_bellymask_WAT, _set_bellymask_WAT)

    def _get_bellymask_ions(self):
        return self._data.get("bellymask_ions", False)
    def _set_bellymask_ions(self, yn):
        self._data["bellymask_ions"] = bool(yn)
    bellymask_ions = property(_get_bellymask_ions, _set_bellymask_ions)


    # other amber parameters
    def set_param(self, param, value):
        self._data.setdefault("cntrl", {})
        param = str(param)
        value = str(value)
        self._data["cntrl"][param] = value

    # ==================================================================
    # properties: filepath
    # ==================================================================

    # for leap
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

    def _get_leap_input_pdb_filepath(self):
        path = self._data.get('leap_input_pdb_filepath', 'leap_input.pdb')
        return path
    leap_input_pdb_filepath = property(_get_leap_input_pdb_filepath)


    # for this task
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
        self._data.setdefault("amber_chain_res_table", {})
        key = "{}:{}".format(amb_chain, amb_res)
        value = "{}:{}".format(orig_chain, orig_res)
        self._data["amber_chain_res_table"][key] = value

    def _get_chainres_match_table(self, amb_chain, amb_res):
        self._data.setdefault("amber_chain_res_table", {})
        key = "{}:{}".format(amb_chain, amb_res)

        value = self._data["amber_chain_res_table"].get(key, None)
        chain = None
        res = None, None
        if value:
            chain, res = value.split(":")
        answer = (chain, res)
        return answer

    # ==================================================================
    # optimization
    # ==================================================================
    def opt(self):
        model = pdfbridge.AtomGroup(self.model)
        self._save_input_pdb(model)

        # remove water?
        model = self._remove_water(model)

        if self.use_belly:
            # belly-mask-atoms are set at top of pdb.
            model = self._reform_model_for_belly(model)

        # self.solvation_model = "cap"

        self._save_leap_input_pdb(model)
        self._prepare_leapin()
        self._do_leap()

        self.set_param("ntx", 1)

        # Every ntpr steps, energy information will be printed
        # in human-readable form to files "mdout" and "mdinfo".
        # default: 50.
        self.set_param("ntpr", 200)

        # Flag for SHAKE to perform bond length constraints.
        # = 1 SHAKE is not performed (default)
        # = 2 bonds involving hydrogen are constrained
        # = 3 allbondsareconstrained(notavailableforparallelorqmmmrunsinsander)
        self.set_param("ntc", 1)

        # Force evaluation.
        # = 1 complete interaction is calculated (default)
        self.set_param("ntf", 1)
        # PBC
        self.set_param("ntb", 0)
        # cutoff range
        self.set_param("cut", 9999.999)
        # Flag for 3D-reference interaction site model
        # idb=1: GB
        self.set_param("igb", 1)

        self._prepare_mdin()

        self._do_sander()
        self._restrt2pdb()


    def _remove_water(self, model):
        return remove_WAT(model)


    def _save_input_pdb(self, model):
        """
        入力ファイルを保存する
        """
        self.cd_workdir("debug input pdb file")

        models = pdfbridge.AtomGroup()
        models.set_group('model_1', model)
        pdb = pdfbridge.Pdb(mode="amber")
        pdb.set_by_atomgroup(models)

        with open(self.input_pdb_filepath, "w") as f:
            f.write(str(pdb))

        self.restore_cwd()

    def _save_leap_input_pdb(self, model):
        """
        leapで処理する前のmodelをpdb形式で保存する
        """
        self.cd_workdir("prepare input pdb file")

        models = pdfbridge.AtomGroup()
        models.set_group('model_1', model)
        pdb = pdfbridge.Pdb(mode="amber")
        pdb.set_by_atomgroup(models)

        with open(self.leap_input_pdb_filepath, "w") as f:
            f.write(str(pdb))

        self.restore_cwd()


    def _reform_model_for_belly(self, model):
        chain_prefix = chr(ord("A") -1)

        if self.bellymask_WAT:
            model_new = pdfbridge.AtomGroup()
            for chain_id, chain in model.groups():
                chain_keep = pdfbridge.AtomGroup()
                chain_keep.name = chain.name
                chain_head = pdfbridge.AtomGroup()
                chain_head.name = chain.name
                for res_id, res in chain.groups():
                    if res.name == "WAT":
                        chain_head.set_group(res_id, res)
                    else:
                        chain_keep.set_group(res_id, res)
                model_new.set_group(chain_id, chain_keep)
                model_new.set_group(chain_prefix + chain_id, chain_head)
            model = self._reorder_model(model_new)

        if self.bellymask_ions:
            ions = ("Na+", "Na",
                    "Cl-", "Cl")
            model_new = pdfbridge.AtomGroup()
            for chain_id, chain in model.groups():
                chain_keep = pdfbridge.AtomGroup()
                chain_keep.name = chain.name
                chain_head = pdfbridge.AtomGroup()
                chain_head.name = chain.name
                for res_id, res in chain.groups():
                    if res.name in ions:
                        chain_head.set_group(res_id, res)
                    else:
                        chain_keep.set_group(res_id, res)
                model_new.set_group(chain_id, chain_keep)
                model_new.set_group(chain_prefix + chain_id, chain_head)
            model = self._reorder_model(model_new)

        return model

    def _reorder_model(self, model):
        input_model = pdfbridge.AtomGroup(model)
        input_model.sort_groups = "nice"

        output_model = pdfbridge.AtomGroup()
        chain_index = 0
        for chain_id, chain in input_model.groups():
            if chain.get_number_of_all_atoms() > 0:
                output_model.set_group(chr(ord("A") + chain_index), chain)
                chain_index += 1

        return output_model


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
        self._save_input_pdb(model)

        # remove water?
        model = self._remove_water(model)

        if self.use_belly:
            # belly-mask-atoms are set at top of pdb.
            model = self._reform_model_for_belly(model)

        # self.solvation_model = "cap"

        self._save_leap_input_pdb(model)
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
        leapin_contents += 'protein = loadPdb {pdb_file}\n'.format(pdb_file=self.leap_input_pdb_filepath)
        leapin_contents += 'proteinBox = copy protein\n'
        leapin_contents += self.__get_leap_ssbond_lines()
        leapin_contents += self.__get_solvation(solute='proteinBox')
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

    def __get_solvation(self, solute):
        """return solvation leap code.
        """
        answer = ''

        if self.solvation_method != None:
            if self.solvation_method == 'cap':
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
                solvent = self.solvation_model
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

        belly_contents = ""
        if self.use_belly:
            pdb = pdfbridge.Pdb()
            pdb.load(self.initial_pdb_filepath)
            initial_model = get_model(pdb.get_atomgroup())
            #print("initial model: #atoms=", initial_model.get_number_of_all_atoms())
            #print("model: ", self.model.get_number_of_all_atoms())
            #print("BM ions", self.bellymask_ions)
            #assert(initial_model.get_number_of_all_atoms() == self.model.get_number_of_all_atoms())

            belly_maskstr = ""
            if self.bellymask_WAT:
                wat_resid = self._get_wat_resid(initial_model)
                if len(belly_maskstr) > 0:
                    belly_maskstr += " | "
                belly_maskstr += self._make_maskstr(resid_areas = wat_resid)
            if self.bellymask_ions:
                ion_resid = self._get_ion_resid(initial_model)
                if len(belly_maskstr) > 0:
                    belly_maskstr += " | "
                belly_maskstr += self._make_maskstr(resid_areas = ion_resid)
            belly_contents = "ibelly=1, bellymask='{belly_maskstr}'"
            belly_contents = belly_contents.format(belly_maskstr=belly_maskstr)


        mdin_contents = """
        #
        &cntrl
          {minimization_contents}
          {cntrl_contents}
          {belly_contents}
        &end
        """

        cntrl_contents = ""
        for param, value in self._data["cntrl"].items():
            cntrl_contents += "{param}={value}, ".format(
                param=param,
                value=value)

        mdin_contents = mdin_contents.format(
            minimization_contents=minimization_contents,
            cntrl_contents=cntrl_contents,
            belly_contents=belly_contents)
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

        # get pdb file from trajectory file using Amber
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

        # match the pdb formed by amber with original one
        orig_model = self._rename_amb2orig(amb_model)
        orig_models = pdfbridge.AtomGroup()
        orig_models.set_group(1, orig_model)
        orig_pdb = pdfbridge.Pdb()
        orig_pdb.set_by_atomgroup(orig_models)
        with open(self.output_pdb_filepath, "w") as f:
            f.write(str(orig_pdb))

        # output for QcModeling
        self.output_model = orig_model

        self.restore_cwd()
        return return_code


    def _rename_amb2orig(self, amb_model):
        """ 対応表(match_table)
        """
        max_chain_id = find_max_chain_id(self.model)
        next_chain_id = chr(ord(max_chain_id) +1)

        # rename model, chain, res_id
        answer = pdfbridge.AtomGroup()
        for amb_chain_id, amb_chain in amb_model.groups():
            for amb_res_id, amb_res in amb_chain.groups():
                orig_chain_id, orig_res_id = self._get_chainres_match_table(amb_chain_id, amb_res_id)
                if (orig_chain_id == None) or (orig_res_id == None):
                    orig_chain_id = next_chain_id
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

        TODO: time-consuming routine!
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
            logger.info('some atoms are deleteed via the leap treatment.: original={} new={}'.format(
                self.model.get_number_of_all_atoms(), amb_model.get_number_of_all_atoms()))

        # 入力モデルとleap後のモデルの座標は変わらないので、
        # 座標をもとに対応表を作成する

        # 座標が等しい原子をmodelから探す
        def find_atom_path(model, atom):
            assert(isinstance(model, pdfbridge.AtomGroup))
            assert(isinstance(atom, pdfbridge.Atom))

            NEAR_DISTANCE = 0.1
            symbol_selector = pdfbridge.Select_Atom(atom.symbol)
            range_selector = pdfbridge.Select_Range(atom.xyz, NEAR_DISTANCE)
            selection = model.select(symbol_selector).select(range_selector)
            path_list = selection.get_path_list()

            answer = None
            if len(path_list) > 0:
                answer = path_list.pop(0)

            return answer

        # make original table
        for chain_id, chain in amb_model.groups():
            for res_id, res in chain.groups():
                for atom_id, atom in res.atoms():
                    amb_path = atom.path
                    orig_path = find_atom_path(self.model, atom)

                    if orig_path == None:
                        pass
                    else:
                        (orig_chain_id, orig_resid, orig_atom_id) = pdfbridge.AtomGroup.divide_path(orig_path)
                        (amb_model_id, amb_chain_id, amb_resid, amb_atom_id) = pdfbridge.AtomGroup.divide_path(amb_path)
                        self._set_chainres_match_table(amb_chain_id, amb_resid, orig_chain_id, orig_resid)

        self._show_amber_chain_res_match_table()

        self.restore_cwd()
        return 0


    def _show_amber_chain_res_match_table(self):
        ''' output matching table
        '''
        logger.debug(">>>> matching table")
        for amb_chain_res, orig_chain_res in self._data["amber_chain_res_table"].items():
            logger.debug("amb: {} <-> orig: {}".format(amb_chain_res, orig_chain_res))


    # ==================================================================
    # for belly method
    # ==================================================================
    def _get_wat_resid(self, model):
        start = -1
        end = -1

        select_WAT = pdfbridge.Select_Name('WAT')
        WATs = model.select(select_WAT)

        return self._get_ambermask_res_list(WATs)


    def _get_ion_resid(self, model):
        start = -1
        end = -1

        areas = []
        select_Na = pdfbridge.Select_Name("Na+")
        select_Cl = pdfbridge.Select_Name("Cl-")
        model_Na = model.select(select_Na)
        model_Cl = model.select(select_Cl)
        model_X = model_Na
        model_X |= model_Cl

        return self._get_ambermask_res_list(model_X)


    def _get_ambermask_res_list(self, select_model):
        select_model.sort_atoms = "nice"
        select_model.sort_groups = "nice"

        start = end = -1
        res_list = []
        for chain_name, chain in select_model.groups():
            for resid, res in chain.groups():
                i = int(resid)
                if start == -1:
                    start = i
                    end = i
                    continue

                if end +1 == i:
                    # continue residue area
                    end = i
                else:
                    res_list.append((start, end))
                    start = i
                    end = i
        if start != -1:
            res_list.append((start, end))

        return res_list


    def _make_maskstr(self, resid_areas):
        assert(isinstance(resid_areas, list))
        maskstr = ""
        for start, end in resid_areas:
            if len(maskstr) == 0:
               maskstr = ":"
            else:
               maskstr += ","

            if start != end:
                maskstr += "{start}-{end}".format(
                    start=start, end=end)
            else:
                maskstr += "{start}".format(
                    start=start)
        return maskstr


if __name__ == '__main__':
    logger.setLevel(logging.DEBUG)
    import doctest
    doctest.testmod()
