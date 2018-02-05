#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path
import shutil
import copy
import math
import yaml
import logging
import msgpack

import proteindf_bridge
from .qccontrol import QcControl

logger = logging.getLogger(__name__)
ANG_PER_AU = 0.5291772

class QcOptRecord(object):
    def __init__(self):

        self.threshold_max_force = 0.00045
        self.threshold_rms_force = 0.00030
        self.threshold_max_displacement = 0.0018
        self.threshold_rms_displacement = 0.0012

        self._steps = []

    def add_step(self, molecule):
        num_of_atoms = molecule.get_number_of_all_atoms()
        coord_table = []
        force_table = []
        self._expand_molecule_data(molecule, coord_table, force_table)
        assert(len(coord_table) == num_of_atoms)
        assert(len(force_table) == num_of_atoms)
        self._steps.append({'molecule': pdfbridge.AtomGroup(molecule),
                            'coord': coord_table,
                            'force': force_table})

    def _expand_molecule_data(self, atomgroup, coord_table, force_table):
        for k, subgroup in atomgroup.groups():
            self._expand_molecule_data(subgroup, coord_table, force_table)
        for k, atom in atomgroup.atoms():
            coord_table.append(atom.xyz)
            force_table.append(atom.force)


    @property
    def step(self):
        step = len(self._steps) -1
        return step

    def _get_stat(self, step):
        assert(0 <= step and step < len(self._steps))
        max_force = 0.0
        rms_force = 0.0
        max_disp = 0.0
        rms_disp = 0.0

        force_table = self._steps[step]['force']
        for force in force_table:
            fx = force.x
            fy = force.y
            fz = force.z
            max_force = max(max_force, abs(fx), abs(fy), abs(fz))
            rms_force += fx*fx + fy*fy + fz*fz
        rms_force = math.sqrt(rms_force / (3.0 * len(force_table)))

        if step > 0:
            coord_table0 = self._steps[step -1]['coord']
            coord_table1 = self._steps[step]['coord']
            assert(len(coord_table0) == len(coord_table1))
            for xyz0, xyz1 in zip(coord_table0, coord_table1):
                x = xyz1.x - xyz0.x
                y = xyz1.y - xyz0.y
                z = xyz1.z - xyz0.z
                max_disp = max(max_disp, abs(x), abs(y), abs(z))
                rms_disp += x*x + y*y + z*z
            rms_disp = math.sqrt(rms_disp / (3.0 * len(coord_table0)))

        return (max_force, rms_force, max_disp, rms_disp)


    def is_converged(self):
        max_force, rms_force, max_disp, rms_disp = self._get_stat(self.step)

        judge_max_force = max_force < self.threshold_max_force
        judge_rms_force = rms_force < self.threshold_rms_force
        judge_max_disp = max_disp < self.threshold_max_displacement
        judge_rms_disp = rms_disp < self.threshold_rms_displacement

        logger.info('MAX FORCE: {} [{}]=> {}'.format(max_force,
                                                     self.threshold_max_force,
                                                     judge_max_force))
        logger.info('RMS FORCE: {} [{}]=> {}'.format(rms_force,
                                                     self.threshold_rms_force,
                                                     judge_rms_force))
        logger.info('MAX DISP.: {} [{}]=> {}'.format(max_disp,
                                                     self.threshold_max_displacement,
                                                     judge_max_disp))
        logger.info('RMS DISP.: {} [{}]=> {}'.format(rms_disp,
                                                     self.threshold_rms_displacement,
                                                     judge_rms_disp))
        judge = judge_max_force and judge_rms_force and judge_max_disp and judge_rms_disp
        return judge


    def get_new_coord(self):
        ans = self.get_SD()
        return ans

    def get_SD(self):
        alpha = 1.0 / (self.step +1)
        molecule = pdfbridge.AtomGroup(self._steps[self.step]['molecule'])

        def update_coord_SD(atomgroup, alpha):
            for k, subgrp in atomgroup.groups():
                update_coord(subgrp, alpha)
            for k, atom in atomgroup.atoms():
                atom.xyz = atom.xyz - ANG_PER_AU * alpha * atom.force
        update_coord_SD(molecule, alpha)
        return molecule


class QcOpt(object):
    def __init__(self,
                 brd_path,
                 template_path =None):

        self._initialize(brd_path, template_path)
        self._opt_record = QcOptRecord()

    @property
    def top_dir(self):
        return self._top_dir

    @property
    def name(self):
        return self._name

    @property
    def step(self):
        return self._step

    @property
    def max_steps(self):
        return self._max_steps

    # --------------------------------------------------------------------------
    def _initialize(self, brd_path, template_path):
        self._top_dir = os.path.abspath(os.getcwd())
        self._name = 'opt'

        self._step = 1
        self._max_steps = 100

        # setup
        self._brd_path = os.path.abspath(brd_path)
        self._molecule = self._load_brd_file(self._brd_path)
        if template_path:
            self._template_path = os.path.abspath(template_path)

    def _load_brd_file(self, path):
        atom_group = None
        with open(path, 'rb') as f:
            raw_dat = msgpack.unpackb(f.read())
            atom_group = pdfbridge.AtomGroup(raw_dat)
        return atom_group

    # --------------------------------------------------------------------------
    def _get_workdir(self, step = None):
        if step == None:
            step = self.step
        wd = os.path.join(self.top_dir,
                          '{name}_{step}'.format(name = self.name,
                                                 step = step))
        return wd

    def _get_brd_path(self, itr = None):
        if itr == None:
            itr = self.step

        if itr == 0:
            path = self._brd_path
        else:
            brd_basename = os.path.basename(self._brd_path)
            path = os.path.join(self._get_workdir(itr),
                                '{}_{}.brd'.format(brd_basename, itr))
        return path

    # --------------------------------------------------------------------------
    def run(self):
        for self._step in range(self._step, self.max_steps):
            self.change_workdir()

            self.calc_single_point()

            if self.is_converged() and self._step > 1:
                logger.info('opt condition is satisfied.')
                break
            else:
                self.update_coord()

        self.restore_workdir()

    def change_workdir(self):
        workdir = self._get_workdir()
        if not os.path.exists(workdir):
            os.mkdir(workdir)
        logger.info('change dir: {}'.format(workdir))
        os.chdir(workdir)

    def restore_workdir(self):
        logger.info('restore dir: {}'.format(self.top_dir))
        os.chdir(self.top_dir)


    def calc_single_point(self):
        logger.info('create QCLO senario')
        yaml_path = self._make_input()

        logger.info('execute QCLO sernario')
        control = QcControl()
        control.run(yaml_path)

        frame = control.get_frame(control.last_frame_name)
        pdfparam = frame.pdfparam
        self._record_step(pdfparam)


    def _make_input(self):
        input_data = {}
        if self._template_path:
            with open(self._template_path, 'r') as file:
                input_data = yaml.load(file)

        # check 'default' task in template
        input_data.setdefault('tasks', [])
        task_default_index = None
        for i, task in enumerate(input_data['tasks']):
            if ('name' in task) and (task['name'] == 'default'):
                task_default_index = i
                break
        if task_default_index == None:
            task_default = {}
            task_default['name'] = 'default'
            input_data['tasks'].append(task_default)
            task_default_index = 0

        # save brd_file in task default
        input_data['tasks'][task_default_index]['brd_file'] = self._get_brd_path(self.step)
        with open(self._get_brd_path(self.step), 'wb') as f:
            f.write(msgpack.packb(self._molecule.get_raw_data()))

        # save YAML file
        yaml_file_path = os.path.join(self._get_workdir(), 'senario.yaml')
        yaml_file = open(yaml_file_path, 'w')
        yaml_file.write(yaml.dump(input_data))
        yaml_file.close

        return yaml_file_path


    def _record_step(self, pdfparam):
        num_of_atoms = pdfparam.num_of_atoms
        logger.info('# of atoms: {}'.format(num_of_atoms))

        forces = []
        for atom_id in range(num_of_atoms):
            force = pdfparam.get_gradient(atom_id)
            forces.append(copy.copy(force))

        num_of_forces = self._insert_forces(self._molecule, forces, 0)
        assert(num_of_forces == num_of_atoms)

        self._opt_record.add_step(self._molecule)


    def _insert_forces(self, atomgroup, force_list, force_list_index):
        for key, group in atomgroup.groups():
            force_list_index = self.insert_forces(group, force_list, force_list_index)
        for key, atom in atomgroup.atoms():
            atom.force = force_list[force_list_index]
            force_list_index += 1
        return force_list_index


    def is_converged(self):
        return self._opt_record.is_converged()

    def update_coord(self):
        self._molecule = self._opt_record.get_new_coord()
