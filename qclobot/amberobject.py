#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil

from .utils import atomgroup2file, file2atomgroup, get_model, check_format_model, find_max_chain_id, remove_WAT
from .process import Process
from .mdobject import MdObject

import proteindf_bridge as bridge

import logging
logger = logging.getLogger(__name__)


class AmberObject(MdObject):
    """Operate Amber task

    : Example

    prepare molecule by AtomGroup(@bridge)
    >>> pdb = bridge.Pdb("./data/sample/1AKG.pdb")
    >>> models = pdb.get_atomgroup()
    >>> model = get_model(models)

    make working directory
    >>> md1 = AmberObject(name='md_test1')

    set model
    >>> md1.model = model

    optimization
    >>> md1.opt()

    """

    def __init__(self, model, work_dir):
        super().__init__(model, work_dir)

    def _initialize(self):
        ''' initialize object
        '''
        super()._initialize()  # called from the parents class

        self._AMBERHOME = os.environ.get('AMBERHOME', '')
        if len(self._AMBERHOME) == 0:
            logger.warning("The environment variable `AMBERHOME` is empty.")

        self._data['leap_sources'] = ['leaprc.protein.ff14SB',
                                      'leaprc.water.tip3p',
                                      'leaprc.gaff']
        self._data['leap_amberparams'] = []

        # default values
        self._is_save_amberparam = True

    # ==================================================================
    # properties
    # ==================================================================
    def _get_is_save_amberparam(self):
        return self._is_save_amberparam

    def _set_is_save_amberparam(self, yn):
        self._is_save_amberparam = bool(yn)
    is_save_amberparam = property(_get_is_save_amberparam,
                                  _set_is_save_amberparam)

    # solvation
    # def _get_is_solvation(self):
    #    return self._data.get('is_solvation', False)
    # def _set_is_solvation(self, yn):
    #    self._data['is_solvation'] = bool(yn)
    # is_solvation = property(_get_is_solvation, _set_is_solvation)

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

    # Generalized Born
    def _get_use_gb(self):
        return self._data.get("use_generalized_born")

    def _set_use_gb(self, yn):
        self._data["use_generalized_born"] = bool(yn)
    use_gb = property(_get_use_gb, _set_use_gb)

    # -------------------------------------------------------------------------
    # restraint
    def _get_use_restraint(self):
        return self._data.get("use_restraint", False)

    def _set_use_restraint(self, yn):
        self._data["use_restraint"] = bool(yn)
    use_restraint = property(_get_use_restraint, _set_use_restraint)

    def _get_restraint_weight(self):
        return self._get_md_param("restraint_wt", 10.0)

    def _set_restraint_weight(self, w):
        self._set_md_param("restraint_wt", w)
    restraint_weight = property(_get_restraint_weight, _set_restraint_weight)

    def _get_restraint_mask(self):
        return self._get_md_param("restraintmask", "")

    def _set_restraint_mask(self, maskstr):
        self._set_md_param("restraintmask", str(maskstr))
    restraint_mask = property(_get_restraint_mask, _set_restraint_mask)

    # -------------------------------------------------------------------------
    # belly type dynamics
    def _get_use_belly(self):
        return self._data.get("use_belly", False)

    def _set_use_belly(self, yn):
        self._data["use_belly"] = bool(yn)
    use_belly = property(_get_use_belly, _set_use_belly)

    def _get_bellymask_H(self):
        return self._data.get("bellymask_H", False)

    def _set_bellymask_H(self, yn):
        self._data["bellymask_H"] = bool(yn)
    bellymask_H = property(_get_bellymask_H, _set_bellymask_H)

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

    def _set_md_param(self, param, value):
        self._data.setdefault("md/cntrl", {})
        self._data["md/cntrl"][param] = value

    def _get_md_param(self, param, default_value=None):
        self._data.setdefault("md/cntrl", {})
        return self._data["md/cntrl"].get(param, default_value)

    def _has_md_param(self, param):
        self._data.setdefault("md/cntrl", {})
        return param in self._data["md/cntrl"]

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

    # def _get_leap_input_pdb_filepath(self):
    #     self._data.setdefault('leap_input_pdb_filepath', 'leap_input.pdb')
    #     return self._data['leap_input_pdb_filepath']
    # leap_input_pdb_filepath = property(_get_leap_input_pdb_filepath)

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
    # complement
    # ==================================================================
    def complement(self):
        """complement protein model
        """
        self._save_input_pdb(self.model)

        self.is_save_amberparam = False
        self._prepare_leapin(input_pdb_filepath=self.input_pdb_filepath,
                             initial_pdb_filepath=self.initial_pdb_filepath)

        answer = self._do_leap()
        self._make_matching_table(self.initial_pdb_filepath)

        # MD doesn't run, so it only runs a copy.
        shutil.copyfile(self.initial_pdb_filepath,
                        self.final_pdb_filepath)

        models = file2atomgroup(self.final_pdb_filepath)
        model = get_model(models)

        renamed_model = self._rename_amb2orig(model)
        self.output_model = renamed_model

        return answer

    # ==================================================================
    # optimization
    # ==================================================================
    def opt(self):
        self._save_input_pdb(self.model)

        # remove water?
        model = self._remove_water(self.model)

        if self.use_belly:
            # belly-mask-atoms are set at top of pdb.
            model = self._reform_model_for_belly(model)

        # self.solvation_model = "cap"

        # self._save_leap_input_pdb(model)
        self._prepare_leapin(input_pdb_filepath=self.input_pdb_filepath,
                             initial_pdb_filepath=self.initial_pdb_filepath)
        self._do_leap()

        self._setup_restraint()

        self._set_md_param("imin", 1)
        self._set_md_param("ntx", 1)
        self._set_md_param("ntpr", 200)
        self._set_md_param("ntc", 1)
        self._set_md_param("ntf", 1)
        self._set_md_param("ntb", 0)
        self._set_md_param("cut", 9999.999)
        # self._set_md_param("igb", 1)

        self._set_md_param("maxcyc", 50000)
        self._set_md_param("ncyc", 0)
        self._set_md_param("ntmin", 3)
        self._set_md_param("drms", 0.001)

        self._prepare_mdin()

        self._do_sander()
        self._restrt2pdb()

        return True

    def _remove_water(self, model):
        return remove_WAT(model)

    def _save_input_pdb(self, model):
        """
        入力ファイルを保存する
        """
        atomgroup2file(model, self.input_pdb_filepath, mode="amber")

    # def _save_leap_input_pdb(self, model):
    #     """
    #     leapで処理する前のmodelをpdb形式で保存する
    #     """
    #     atomgroup2file(model, self.leap_input_pdb_filepath)
        # self.cd_workdir("prepare input pdb file")

        # models = bridge.AtomGroup()
        # models.set_group('model_1', model)
        # pdb = bridge.Pdb(mode="amber")
        # pdb.set_by_atomgroup(models)

        # with open(self.leap_input_pdb_filepath, "w") as f:
        #     f.write(str(pdb))

        # self.restore_cwd()

    def _setup_restraint(self):
        if self.use_restraint:
            self._set_md_param("ntr", 1)

        if self.use_belly:
            bellymask_str = ""
            if self.bellymask_WAT:
                # wat_resid = self._get_wat_resid(initial_model)
                if len(bellymask_str) > 0:
                    bellymask_str += " | "
                # belly_maskstr += self._make_maskstr(resid_areas = wat_resid)
                bellymask_str += ":WAT"
            if self.bellymask_ions:
                # ion_resid = self._get_ion_resid(initial_model)
                if len(bellymask_str) > 0:
                    bellymask_str += " | "
                # belly_maskstr += self._make_maskstr(resid_areas = ion_resid)
                bellymask_str += ":Na+,Cl-"
            if self.bellymask_H:
                if len(bellymask_str) > 0:
                    bellymask_str += " | "
                bellymask_str += "@H="

            self._set_md_param("ibelly", 1)
            self._set_md_param("bellymask", bellymask_str)

    def _reform_model_for_belly(self, model):
        chain_prefix = chr(ord("A") - 1)

        if self.bellymask_WAT:
            model_new = bridge.AtomGroup()
            for chain_id, chain in model.groups():
                chain_keep = bridge.AtomGroup()
                chain_keep.name = chain.name
                chain_head = bridge.AtomGroup()
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
            model_new = bridge.AtomGroup()
            for chain_id, chain in model.groups():
                chain_keep = bridge.AtomGroup()
                chain_keep.name = chain.name
                chain_head = bridge.AtomGroup()
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
        input_model = bridge.AtomGroup(model)
        input_model.sort_groups = "nice"

        output_model = bridge.AtomGroup()
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
           steps=1, dt=0.002):
        """ compute molecular dynamics

        steps: Number of MD-steps to be performed.
        dt:    The time step (psec).
        """
        model = bridge.AtomGroup(self.model)
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

        self._set_md_param("imin", 0)
        self._set_md_param("irest", 0)

        self._set_md_param("ntpr", 5000)
        self._set_md_param("ntwx", 500)
        self._set_md_param("ntwe", 500)

        self._set_md_param("nstlim", steps)  # 500000
        self._set_md_param("dt", dt)

        self._set_md_param("ntt", 1)
        self._set_md_param("temp0", 300.0)
        self._set_md_param("tempi", 0.0)
        self._set_md_param("tautp", 0.1)

        # SHAKE bond length constraints
        self._set_md_param("ntc", 2)

        self._set_md_param("ntf", 2)
        self._set_md_param("ntb", 0)
        self._set_md_param("cut", 10.0)

        self._prepare_mdin(opt=False)

        self._do_sander()
        self._restrt2pdb()

    # ==================================================================
    # leap setting
    # ==================================================================

    def _prepare_leapin(self,
                        input_pdb_filepath=None,
                        initial_pdb_filepath=None):
        """make input file for xleap command
        """
        if input_pdb_filepath is None:
            input_pdb_filepath = self.input_pdb_filepath
        if initial_pdb_filepath is None:
            initial_pdb_filepath = self.initial_pdb_filepath

        # self.cd_workdir("prepare leapin")

        leapin_contents = ''
        leapin_contents += 'logFile {logfile_filepath}\n'.format(
            logfile_filepath=self.leap_logfile_filepath)
        leapin_contents += self.__get_leap_source_lines()
        leapin_contents += self.__get_leap_amberparams_lines()
        leapin_contents += 'protein = loadPdb {pdb_file}\n'.format(
            pdb_file=self.input_pdb_filepath)
        leapin_contents += 'proteinBox = copy protein\n'
        leapin_contents += self.__get_leap_ssbond_lines()
        leapin_contents += self.__get_solvation(solute='proteinBox')

        if self.is_save_amberparam:
            leapin_contents += 'saveAmberParm proteinBox {prmtop} {inpcrd}\n'.format(prmtop=self.prmtop_filepath,
                                                                                     inpcrd=self.inpcrd_filepath)
        leapin_contents += 'savePdb proteinBox {pdb_filepath}\n'.format(
            pdb_filepath=initial_pdb_filepath)
        leapin_contents += 'quit\n'

        logger.debug('save leap inputfile: {}'.format(self.leapin_filepath))
        with open(self.leapin_filepath, 'w') as f:
            f.write(leapin_contents)

        # self.restore_cwd()

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

        if self.solvation_method is not None:
            if self.solvation_method == 'cap':
                center = self.model.center()
                cap_center = "{{ {x:.3f} {y:.3f} {z:.3f} }}".format(x=center.x,
                                                                    y=center.y,
                                                                    z=center.z)
                (box_min, box_max) = self.model.box()
                distance = 0.0
                distance = max(distance, center.distance_from(
                    bridge.Position(box_min.x, box_min.y, box_min.z)))
                distance = max(distance, center.distance_from(
                    bridge.Position(box_max.x, box_min.y, box_min.z)))
                distance = max(distance, center.distance_from(
                    bridge.Position(box_max.x, box_max.y, box_min.z)))
                distance = max(distance, center.distance_from(
                    bridge.Position(box_max.x, box_min.y, box_max.z)))
                distance = max(distance, center.distance_from(
                    bridge.Position(box_max.x, box_max.y, box_max.z)))
                distance = max(distance, center.distance_from(
                    bridge.Position(box_min.x, box_max.y, box_min.z)))
                distance = max(distance, center.distance_from(
                    bridge.Position(box_min.x, box_max.y, box_max.z)))
                distance = max(distance, center.distance_from(
                    bridge.Position(box_min.x, box_min.y, box_max.z)))
                cap_radius = max(distance + 10.0, 30.0)
                solvent = self.solvation_model
                answer += "solvateCap {solute} {solvent} {position} {radius} {closeness}\n".format(
                    solute=solute,
                    solvent=solvent,
                    position=cap_center,
                    radius=cap_radius,
                    closeness="")
            else:
                logger.warning('solvation model is not understood.: {}'.format(
                    self.solvation_model))

        return answer

    # ==================================================================
    # do leap
    # ==================================================================
    def _do_leap(self):
        # elf.cd_workdir("do leap")

        p = Process()
        leap_cmd = os.path.join(self._AMBERHOME, 'bin', 'tleap')
        cmd = "{leap_cmd} -s -f {leapin}".format(leap_cmd=leap_cmd,
                                                 leapin=self.leapin_filepath)

        logger.info("run command: {}".format(cmd))
        p.cmd(cmd)
        return_code = p.commit_blocked('stdout.txt', 'stderr.txt',
                                       stdout_through=False,
                                       stderr_through=False)

        # self.restore_cwd()
        return return_code

    # ==================================================================
    # prepare mdin
    # ==================================================================
    def _prepare_mdin(self):

        def get_line(cmd, remark, indent=2):
            line = ""
            line += " " * indent
            line += cmd + ","
            line += " " * (20 - len(cmd))
            line += "! " + remark + "\n"

            return line

        options = ""
        # 19.6. General minimization and dynamics parameters
        # 19.6.1. General flags describing the calculation
        if self._has_md_param("imin"):
            options += get_line("imin={}".format(self._get_md_param("imin")),
                                "Flag to run minimization")

        # 19.6.2. Nature and format of the input
        if self._has_md_param("ntx"):
            options += get_line("ntx={}".format(self._get_md_param("ntx")),
                                "Option to read the coordinates from the \"inpcrd\" file")

        # 19.6.3. Nature and format of the output
        if self._has_md_param("ntpr"):
            options += get_line("ntpr={}".format(self._get_md_param("ntpr")),
                                "Every ntpr steps, energy information will be printed in human-readable form to files \"mdout\" and \"mdinfo\"")

        # 19.6.4. Frozen or restrained atoms
        if self._has_md_param("ibelly"):
            options += get_line("ibelly={}".format(self._get_md_param("ibelly")),
                                "Flag for belly type dynamics")
            options += get_line("bellymask='{}'".format(self._get_md_param("bellymask")), "")

        if self._has_md_param("ntr"):
            options += get_line("ntr=1",
                                "Flag for restraining specified atoms in Cartesian space using a harmonic potential, if ntr > 0")
            options += get_line("restraint_wt={}".format(self.restraint_weight),
                                "The weight (in kcal/mol−Å2) for the positional restraints. ")
            options += get_line("restraintmask='{}'".format(self._get_md_param("restraintmask")), "")

        # 19.6.5. Energy minimization
        if self._has_md_param("maxcyc"):
            options += get_line("maxcyc={}".format(self._get_md_param("maxcyc")),
                                "The maximum number of cycles of minimization")
        if self._has_md_param("ncyc"):
            options += get_line("ncyc={}".format(self._get_md_param("ncyc")),
                                "If NTMIN is 1 then the method of minimization will be switched from steepest descent to conjugate gradient after NCYC cycles")
        if self._has_md_param("ntmin"):
            options += get_line("ntmin={}".format(self._get_md_param("ntmin")),
                                "Flag for the minimization")
        if self._has_md_param("drms"):
            options += get_line("drms={}".format(self._get_md_param("drms")),
                                "The convergence criterion for the energy gradient")

        # 19.6.9. SHAKE bond length constraints
        if self._has_md_param("ntc"):
            options += get_line("ntc={}".format(self._get_md_param("ntc")),
                                "Flag for SHAKE to perform bond length constraints")

        # 19.7. Potential function parameters
        # 19.7.1. Generic parameters
        if self._has_md_param("ntf"):
            options += get_line("ntf={}".format(self._get_md_param("ntf")), "Force evaluation")
        if self._has_md_param("ntb"):
            options += get_line("ntb={}".format(self._get_md_param("ntb")),
                                "This variable controls whether or not periodic boundaries are imposed on the system during the calculation of non-bonded interactions")
        if self._has_md_param("cut"):
            options += get_line("cut={}".format(self._get_md_param("cut")),
                                "This is used to specify the nonbonded cutoff, in Angstroms")

        with open(self.mdin_filepath, 'w') as f:
            f.write("#\n")
            f.write("&cntrl\n")
            f.write(options)
            f.write("&end\n")

    # ==================================================================
    # do sander
    # ==================================================================
    def _do_sander(self):
        # self.cd_workdir("do sander")

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

        if self.use_restraint:
            cmd += " -ref {refc}".format(refc=self.inpcrd_filepath)

        logger.info("run command: {}".format(cmd))
        p.cmd(cmd)
        return_code = p.commit()

        # self.restore_cwd()
        logger.info("sander done.")
        return return_code

    # ==================================================================
    # get pdb
    # ==================================================================
    def _restrt2pdb(self):

        self._make_matching_table()

        # self.cd_workdir("restrt2pdb")

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
        logger.info("save amber pdb file: {}".format(self.final_pdb_filepath))

        amb_pdb = bridge.Pdb(self.final_pdb_filepath)
        amb_models = amb_pdb.get_atomgroup()
        amb_model = get_model(amb_models)

        # match the pdb formed by amber with original one
        orig_model = self._rename_amb2orig(amb_model)
        orig_models = bridge.AtomGroup()
        orig_models.set_group(1, orig_model)
        orig_pdb = bridge.Pdb()
        orig_pdb.set_by_atomgroup(orig_models)
        logger.info("save reordered pdb file: {}".format(
            self.output_pdb_filepath))
        with open(self.output_pdb_filepath, "w") as f:
            f.write(str(orig_pdb))

        # output for QcModeling
        self.output_model = orig_model

        # self.restore_cwd()
        return return_code

    def _rename_amb2orig(self, amb_model):
        """ 対応表(match_table)
        """
        max_chain_id = find_max_chain_id(self.model)
        next_chain_id = chr(ord(max_chain_id) + 1)

        # rename model, chain, res_id
        answer = bridge.AtomGroup()
        for amb_chain_id, amb_chain in amb_model.groups():
            for amb_res_id, amb_res in amb_chain.groups():
                orig_chain_id, orig_res_id = self._get_chainres_match_table(
                    amb_chain_id, amb_res_id)
                logger.debug("amb: {amb_chain_id}/{amb_res_id} -> new: {orig_chain_id}/{orig_res_id}".format(
                    amb_chain_id=amb_chain_id, amb_res_id=amb_res_id, orig_chain_id=orig_chain_id, orig_res_id=orig_res_id))

                if (orig_chain_id is None) or (orig_res_id is None):
                    # logger.warning("Not found amber ID: {}/{}".format(amb_chain_id, amb_res_id))
                    orig_chain_id = next_chain_id
                    orig_res_id = amb_res_id

                if answer.has_groupkey(orig_chain_id) is not True:
                    orig_chain = bridge.AtomGroup(name=amb_chain.name)
                    answer.set_group(orig_chain_id, orig_chain)
                answer[orig_chain_id].set_group(orig_res_id, amb_res)

        return answer

    # ==================================================================
    # make matching table
    # ==================================================================
    def _make_matching_table(self, initial_pdb_filepath=None):
        """
        make matching table between original and leap model.

        入力したモデルとleapで処理したモデルで残基番号などが変化するので、
        対応表を作成する

        TODO: time-consuming routine!
        """
        if initial_pdb_filepath is None:
            initial_pdb_filepath = self.initial_pdb_filepath

        logger.info("begin: make matching table")
        # self.cd_workdir("match models")

        # make model_amber
        amb_pdb = bridge.Pdb(initial_pdb_filepath)
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
        def find_atom_path(model, atom, symbol_cache):
            assert(isinstance(model, bridge.AtomGroup))
            assert(isinstance(atom, bridge.Atom))
            assert(isinstance(symbol_cache, dict))

            NEAR_DISTANCE = 0.1
            if atom.symbol not in symbol_cache:
                symbol_selector = bridge.Select_Symbol(atom.symbol)
                symbol_cache[atom.symbol] = model.select(symbol_selector)

            range_selector = bridge.Select_Range(atom.xyz, NEAR_DISTANCE)

            # selection1 = model.select(symbol_selector)
            # print(atom.symbol, str(selection1))

            # selection = model.select(symbol_selector).select(range_selector)
            selection = symbol_cache[atom.symbol].select(range_selector)
            path_list = selection.get_path_list()

            answer = None
            if len(path_list) > 0:
                answer = path_list.pop(0)

            return answer

        # make original table
        # print(self.model)
        self._data.setdefault("amber_chain_res_table", {})
        cache = {}
        for chain_id, chain in amb_model.groups():
            for res_id, res in chain.groups():
                for atom_id, atom in res.atoms():
                    amb_path = atom.path
                    orig_path = find_atom_path(self.model, atom, cache)

                    # print("{} <-> {}".format(amb_path, orig_path))
                    if orig_path is None:
                        pass
                    else:
                        (orig_chain_id, orig_resid,
                         orig_atom_id) = bridge.AtomGroup.divide_path(orig_path)
                        (amb_model_id, amb_chain_id, amb_resid,
                         amb_atom_id) = bridge.AtomGroup.divide_path(amb_path)
                        self._set_chainres_match_table(
                            amb_chain_id, amb_resid, orig_chain_id, orig_resid)

        self._show_amber_chain_res_match_table()

        # self.restore_cwd()
        logger.info("end: make matching table")
        return 0

    def _set_chainres_match_table(self,
                                  amb_chain, amb_res,
                                  orig_chain, orig_res):
        key = "{}:{}".format(amb_chain, amb_res)
        value = "{}:{}".format(orig_chain, orig_res)
        self._data["amber_chain_res_table"][key] = value

        resid_only_key = "*:{}".format(amb_res)
        self._data["amber_chain_res_table"].setdefault(resid_only_key, value)

    def _get_chainres_match_table(self, amb_chain, amb_res):
        self._data.setdefault("amber_chain_res_table", {})
        key = "{}:{}".format(amb_chain, amb_res)

        chain = None
        res = None
        value = self._data["amber_chain_res_table"].get(key, None)
        if value is None:
            resid_only_key = "*:{}".format(amb_res)
            value = self._data["amber_chain_res_table"].get(
                resid_only_key, None)
        if value:
            chain, res = value.split(":")
        answer = (chain, res)
        return answer

    def _show_amber_chain_res_match_table(self):
        ''' output matching table
        '''
        logger.debug(">>>> matching table")
        if len(self._data["amber_chain_res_table"]) == 0:
            logger.warning("empty amber chain:res matching table")
        for amb_chain_res, orig_chain_res in self._data["amber_chain_res_table"].items():
            logger.debug(
                "amb: {} <-> orig: {}".format(amb_chain_res, orig_chain_res))

    # ==================================================================
    # for belly method
    # ==================================================================

    def _get_wat_resid(self, model):
        # start = -1
        # end = -1

        select_WAT = bridge.Select_Name('WAT')
        WATs = model.select(select_WAT)

        return self._get_ambermask_res_list(WATs)

    def _get_ion_resid(self, model):
        # start = -1
        # end = -1

        # areas = []
        select_Na = bridge.Select_Name("Na+")
        select_Cl = bridge.Select_Name("Cl-")
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

                if end + 1 == i:
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
