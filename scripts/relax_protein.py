#!/usrbin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 The ProteinDF development team.
# see also AUTHORS and README if provided.
#
# This file is a part of the ProteinDF software package.
#
# The ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import os
import string
import argparse
import subprocess
import shlex
import re

import proteindf_bridge as bridge
import proteindf_tools as pdf
import qclobot as qclo

import logging


def main():
    """
    """
    # parse args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('FILE',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('step',
                        nargs=1,
                        type=int
                        )
    parser.add_argument('-l', '--logfile',
                        nargs=1,
                        action='store',
                        help='logconfig file')
    parser.add_argument('-L', '--logconfig',
                        nargs=1,
                        action='store',
                        help='logconfig file')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    # setting
    mpac_file_path = args.FILE[0]
    step = args.step[0]
    verbose = args.verbose
    is_debug = args.debug

    if args.logconfig:
        # logging module setup by config-file
        logconfig_path = args.logconfig[0]
        logging.config.fileConfig(logconfig_path)
    else:
        # logging module setup
        logfile_path = ''
        if args.logfile:
            logfile_path = args.logfile[0]
        setup_logging(logfile_path, is_debug)

    # reading
    if verbose:
        print("reading: %s\n" % (mpac_file_path))

    # prepare atomgroup
    models = bridge.load_atomgroup(mpac_file_path)
    # print(models)

    #
    r = Relax(models)
    if step == 1:
        r.step1()
    elif step == 2:
        r.step2()
    elif step == 3:
        r.step3()
    else:
        raise

    print('finish')
    exit


def setup_logging(logfile_path='', is_debug=False):
    if len(logfile_path) == 0:
        logfile_path = 'relax.log'

    logging_level = logging.INFO
    format_str = '%(asctime)s [%(levelname)s] %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'
    if is_debug:
        logging_level = logging.DEBUG
        format_str = '%(asctime)s [%(levelname)s] [%(name)s] %(message)s'

    logging.basicConfig(
        filename=logfile_path,
        level=logging_level,
        format=format_str,
        datefmt=date_format
    )

    formatter = logging.Formatter(format_str, date_format)

    console = logging.StreamHandler()
    console.setLevel(logging.WARNING)
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)


class Relax(object):
    def __init__(self, atomgroup):
        self._logger = logging.getLogger(__name__)
        # self._logger.addHandler(logging.NullHandler())
        # self._logger.setLevel(logging.INFO)

        self._chain_residues = None
        self._registered_group = None  # leapで処理できる残基名リスト
        self._res_set = None  # 計算対象の残基名リスト
        self._leapin_premd_filepath = 'leap_pre.in'
        self._leapin_md1_filepath = 'leap_md1.in'
        self._leapin_md2_filepath = 'leap_md2.in'
        self._leapin_md3_filepath = 'leap_md3.in'
        self._leapin_ssbond_filepath = 'leap_ssbond.in'
        self._antechamber_cmd = 'antechamber'
        self._parmchk_cmd = 'parmchk'

        self._leaprc = [
            "leaprc.f14SB",
            "leaprc.protein.ff14SB",
            "leaprc.water.tip3p",
            "leaprc.gaff"
        ]

        self._ligands = []

        # for Amber modpdb
        #biopdb = bridge.Pdb(mode='amber')
        #atomgroup = biopdb.get_modpdb_atomgroup(atomgroup)

        # check models
        self._num_of_models = atomgroup.get_number_of_groups()
        if self._num_of_models > 1:
            self._logger.warn(
                '# of models(={}) > 1'.format(self._num_of_models))
        (self._model_name, self._model) = list(atomgroup.groups())[0]
        self._logger.info("model: {}".format(self._model_name))

        # check model
        self._num_of_chains = self._model.get_number_of_groups()
        self._logger.info('# of chains: {}'.format(self._num_of_chains))

        self._prepare()

    def _prepare(self):
        self._trans_HIS()

        self._logger.info('check groups...')
        self._make_group_DB()
        self._make_res_DB()

        for res in self._res_set:
            if res not in self._registered_group:
                self._logger.info(
                    '{} is not found in registered group.'.format(res))
                self._ligands.append(res)
        self._logger.info('check groups. done.')

        self._logger.debug(str(self._ligands))
        for ligand in self._ligands:
            frcmod_path = '{ligand}.frcmod'.format(ligand=ligand)
            prep_path = '{ligand}.prep'.format(ligand=ligand)
            if ((os.path.exists(frcmod_path) == True) and
                    (os.path.exists(prep_path) == True)):
                self._logger.info(
                    'already created parameters: {ligand}'.format(ligand=ligand))
            else:
                self._make_prep(ligand)

    def step1(self):
        self._prepare_leapin_step1()
        self._do_leap(self._leapin_md1_filepath)
        self._prepare_mdin_step1()

    def step2(self):
        self._prepare_leapin_step2()
        self._do_leap(self._leapin_md2_filepath)
        self._prepare_mdin_step2()

    def step3(self):
        self._prepare_leapin_step3()
        self._do_leap(self._leapin_md3_filepath)
        self._prepare_mdin_step3()

    def _prepare_leapin_step1(self):
        md1_input_pdb_filepath = 'md1_input.pdb'

        leapin_templ = """
        ${source_leaprc}

        # for user-defined parameters
        ${load_ambparam_str}

        protein = loadPdb ${pdb_file}

        ${ssbond_str}

        desc protein
        charge protein
        check protein

        saveAmberParm protein md1.prmtop md1.inpcrd
        savepdb protein md1_before.pdb

        quit
        """
        leapin_templ = leapin_templ.lstrip('\n')
        leapin_templ = bridge.Utils.unindent_block(leapin_templ)

        load_ambparam_str = self._get_load_ambparam_str()
        ssbond_str = self._get_SS_bond_cmd()  # <- update self._model object
        leapin_contents = string.Template(leapin_templ).substitute({
            'source_leaprc': self._get_leaprc_source_lines(),
            'load_ambparam_str': load_ambparam_str,
            'pdb_file': md1_input_pdb_filepath,
            'ssbond_str': ssbond_str
        })

        # output
        self._logger.info(
            'save pdb file for input : {}'.format(md1_input_pdb_filepath))
        self._save_pdb(md1_input_pdb_filepath)

        self._logger.info('save leap.in file: {}'.format(
            self._leapin_md1_filepath))
        fout = open(self._leapin_md1_filepath, 'w')
        fout.write(leapin_contents)
        fout.close()

    def _prepare_leapin_step2(self):
        self._neutralize()
        self._check_ionpairs()

        self._save_pdb('md2_neutralize.pdb')
        self._model = self._reorder_ions_for_amber(self._model)
        self._logger.debug(str(self._model))

        md2_input_pdb_filepath = 'md2_input.pdb'
        self._save_pdb(md2_input_pdb_filepath)

        center = self._model.center()
        solcap_center = "{{ {x:.3f} {y:.3f} {z:.3f} }}".format(x=center.x,
                                                               y=center.y,
                                                               z=center.z)

        (box_min, box_max) = self._model.box()
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
        solcap_closeness = max(distance + 10.0, 30.0)

        leapin_templ = """
        ${source_leaprc}

        # for Na+ and Cl-
        loadamberparams frcmod.ionsjc_tip3p

        # for user-defined parameters
        ${load_ambparam_str}

        protein = loadPdb ${pdb_file}

        ${ssbond_str}

        desc protein
        charge protein
        check protein

        solvateCap protein TIP3PBOX ${solcap_center} ${solcap_closeness}

        saveAmberParm protein md2.prmtop md2.inpcrd
        savepdb protein md2_before.pdb

        quit
        """
        leapin_templ = leapin_templ.lstrip('\n')
        leapin_templ = bridge.Utils.unindent_block(leapin_templ)

        load_ambparam_str = self._get_load_ambparam_str()
        ssbond_str = self._get_SS_bond_cmd()
        leapin_contents = string.Template(leapin_templ).substitute({
            'source_leaprc': self._get_leaprc_source_lines(),
            'load_ambparam_str': load_ambparam_str,
            'pdb_file': md2_input_pdb_filepath,
            'solcap_center': solcap_center,
            'solcap_closeness': solcap_closeness,
            'ssbond_str': ssbond_str
        })

        #leapin_templ = leapin_templ.lstrip('\n')
        #leapin_templ = bridge.Utils.unindent_block(leapin_templ)

        # output
        self._logger.info('save leap.in file: {}'.format(
            self._leapin_md2_filepath))
        fout = open(self._leapin_md2_filepath, 'w')
        fout.write(leapin_contents)
        fout.close()

    def _prepare_leapin_step3(self):
        md3_input_pdb_filepath = 'md3_input.pdb'
        self._save_pdb(md3_input_pdb_filepath)

        leapin_templ = """
        ${source_leaprc}

        # for Na+ and Cl-
        loadamberparams frcmod.ionsjc_tip3p

        # for user-defined parameters
        ${load_ambparam_str}

        protein = loadPdb ${pdb_file}

        ${ssbond_str}

        desc protein
        charge protein
        check protein

        saveAmberParm protein md3.prmtop md3.inpcrd
        savepdb protein md3_before.pdb

        quit
        """
        leapin_templ = leapin_templ.lstrip('\n')
        leapin_templ = bridge.Utils.unindent_block(leapin_templ)

        load_ambparam_str = self._get_load_ambparam_str()
        ssbond_str = self._get_SS_bond_cmd()
        leapin_contents = string.Template(leapin_templ).substitute({
            'source_leaprc': self._get_leaprc_source_lines(),
            'load_ambparam_str': load_ambparam_str,
            'pdb_file': md3_input_pdb_filepath,
            'ssbond_str': ssbond_str
        })

        # output
        self._logger.info('save leap.in file: {}'.format(
            self._leapin_md3_filepath))
        fout = open(self._leapin_md3_filepath, 'w')
        fout.write(leapin_contents)
        fout.close()

    def _save_pdb(self, filepath):
        models = bridge.AtomGroup()
        models[self._model_name] = self._model

        # prepare BrPdb object
        pdb_obj = bridge.Pdb(mode='amber')
        pdb_obj.set_by_atomgroup(models)

        # output PDB
        self._logger.info('output PDB file: {}'.format(filepath))
        fout = open(filepath, "w")
        contents = str(pdb_obj)
        fout.write(contents)
        fout.close()

    def _get_leaprc_source_lines(self):
        """
        """
        answer = ""
        if isinstance(self._leaprc, list):
            for obj in self._leaprc:
                answer += "source {}\n".format(obj)
        else:
            answer = "source {}\n".format(self._leaprc)

        return answer

    def _get_load_ambparam_str(self):
        """
        """
        retval = ''

        for lig in self._ligands:
            retval += 'loadAmberParams {lig}.frcmod\n'.format(lig=lig)
            retval += 'loadAmberPrep {lig}.prep\n'.format(lig=lig)

        return retval

    def _get_SS_bond_cmd(self):
        leap_ssbond_cmd = ''

        if os.path.exists(self._leapin_ssbond_filepath):
            f = open(self._leapin_ssbond_filepath, 'r')
            leap_ssbond_cmd = f.read()
            f.close()
        else:
            leap_ssbond_cmd = ''
            chain_residues = self._get_chain_residues_array()

            ss_bond_residues_pairset = set()
            bond_list = self._model.get_bond_list()
            for bond in bond_list:
                (model_name1, chain_name1, resid1,
                 atom_name1) = self._divide_path(bond[0])
                (model_name2, chain_name2, resid2,
                 atom_name2) = self._divide_path(bond[1])
                self._model[chain_name1][resid1].name = "CYX"
                self._model[chain_name2][resid2].name = "CYX"

                resid1 = self._get_serial_resid(
                    chain_residues, chain_name1, resid1)
                resid2 = self._get_serial_resid(
                    chain_residues, chain_name2, resid2)
                if resid1 > resid2:
                    resid1, resid2 = resid2, resid1  # swap
                ss_bond_residues_pairset.add((resid1, resid2))

            for (resid1, resid2) in ss_bond_residues_pairset:
                leap_ssbond_cmd += 'bond protein.{resid1}.SG protein.{resid2}.SG\n'.format(resid1=resid1,
                                                                                           resid2=resid2)
            f = open(self._leapin_ssbond_filepath, 'w')
            f.write(leap_ssbond_cmd)
            f.close()

        return leap_ssbond_cmd

    def _get_chain_residues_array(self):
        if self._chain_residues == None:
            chain_residues = []
            for chain_name, chain in self._model.groups():
                chain_residues.append(
                    (chain_name, int(chain.get_number_of_groups())))
                self._logger.info('# of residues in chain {chain_name}: {num_res}'.format(chain_name=chain_name,
                                                                                          num_res=chain.get_number_of_groups()))
            self._chain_residues = chain_residues

        return self._chain_residues

    def _divide_path(self, path):
        """
        atomgroupのパスを構成要素ごとに分割する
        """
        items = path.rsplit('/')
        answer = []
        for item in items:
            if len(item) > 0:
                answer.append(item)
        return answer

    def _get_serial_resid(self, chain_residues, chain_name, resid):
        """
        鎖を無視した通し番号を返す
        """
        resid = int(resid)

        index = 0
        for c, n in chain_residues:
            if c == chain_name:
                break
            index += int(n)
        index += resid

        return index

    def _trans_HIS(self):
        self._logger.info('check HIS...')
        for chain_name, chain in self._model.groups():
            for resid, res in chain.groups():
                resname = res.name
                if resname == 'HIS':
                    has_HE2 = False
                    has_HD2 = False
                    if res.has_atomname('HE2') == True:
                        self._logger.info('HE2 found: {}'.format(res.path))
                        has_HE2 = True
                        res.name = 'HIE'
                    elif res.has_atomname('HD2') == True:
                        self._logger.info('HD2 found: {}'.format(res.path))
                        has_HD2 = True
                        res.name = 'HID'

                    if has_HE2 and has_HD2:
                        res.name = 'HIP'
                    elif has_HE2:
                        res.name = 'HIE'
                    elif has_HD2:
                        res.name = 'HID'
                    else:
                        self._logger.warn(
                            'cannot assign HIS: {}/{}'.format(chain_name, resid))
                        for atmkey, atm in res.atoms():
                            self._logger.warn(
                                'atom name: "{}"'.format(atm.name))
        self._logger.info('check HIS. done.')

    def _neutralize(self):
        ip = bridge.IonPair(self._model)
        ionpairs = ip.get_ion_pairs()

        # 処理しやすいように並べ替え
        exempt_list = []
        for (anion_path, cation_path, anion_type, cation_type) in ionpairs:
            (anion_chain_name, anion_res_name) = self._divide_path(anion_path)
            (cation_chain_name, cation_res_name) = self._divide_path(cation_path)
            exempt_list.append((anion_chain_name, anion_res_name, anion_type))
            exempt_list.append(
                (cation_chain_name, cation_res_name, cation_type))

        modeling = bridge.Modeling()
        for chain_name, chain in self._model.groups():
            for resid, res in chain.groups():
                resname = res.name

                if res.has_atom('H3'):
                    if (chain_name, resid, 'NTM') not in exempt_list:
                        # N-term
                        ag = modeling.neutralize_Nterm(res)
                        self._logger.info("add ion for N-term: {}".format(ag))
                        # self._add_ions_to_new_chain(ag)
                        self._add_ions(res, ag)
                    else:
                        self._logger.info(
                            'exempt adding ion: {}/{} Nterm'.format(chain_name, resname))
                if res.has_atom('OXT'):
                    if (chain_name, resid, 'CTM') not in exempt_list:
                        # C-term
                        ag = modeling.neutralize_Cterm(res)
                        self._logger.info("add ion for C-term: {}".format(ag))
                        # self._add_ions_to_new_chain(ag)
                        self._add_ions(res, ag)
                    else:
                        self._logger.info(
                            'exempt adding ion: {}/{} Cterm'.format(chain_name, resname))

                if resname == 'GLU':
                    if (chain_name, resid, 'GLU') not in exempt_list:
                        ag = modeling.neutralize_GLU(res)
                        self._logger.info(
                            "add ion for GLU({}): {}".format(resid, ag))
                        # self._add_ions_to_new_chain(ag)
                        self._add_ions(res, ag)
                    else:
                        self._logger.info(
                            'exempt adding ion: {}/{} GLU'.format(chain_name, resname))
                elif resname == 'ASP':
                    if (chain_name, resid, 'ASP') not in exempt_list:
                        ag = modeling.neutralize_ASP(res)
                        self._logger.info(
                            "add ion for ASP({}): {}".format(resid, ag))
                        # self._add_ions_to_new_chain(ag)
                        self._add_ions(res, ag)
                    else:
                        self._logger.info(
                            'exempt adding ion: {}/{} ASP'.format(chain_name, resname))
                elif resname == 'LYS':
                    if (chain_name, resid, 'LYS') not in exempt_list:
                        ag = modeling.neutralize_LYS(res)
                        self._logger.info(
                            "add ion for LYS({}): {}".format(resid, ag))
                        # self._add_ions_to_new_chain(ag)
                        self._add_ions(res, ag)
                    else:
                        self._logger.info(
                            'exempt adding ion: {}/{} LYS'.format(chain_name, resname))
                elif resname == 'ARG':
                    if (((chain_name, resid, 'ARG') not in exempt_list) and
                        ((chain_name, resid, 'ARG1') not in exempt_list) and
                            ((chain_name, resid, 'ARG2') not in exempt_list)):
                        ag = modeling.neutralize_ARG(res)
                        self._logger.info(
                            "add ion for ARG({}): {}".format(resid, ag))
                        # self._add_ions_to_new_chain(ag)
                        self._add_ions(res, ag)
                    else:
                        self._logger.info(
                            'exempt adding ion: {}/{} ARG'.format(chain_name, resname))

                elif resname == 'FAD':
                    ag = modeling.neutralize_FAD(res)
                    self._logger.info(
                        "add ion for FAD({}): {}".format(resid, ag))
                    # self._add_ions_to_new_chain(ag)
                    self._add_ions(res, ag)

    # def _add_ions_to_new_chain(self, atomgroup):
    #    # for Amber PDB format
    #    chain = bridge.AtomGroup()
    #
    #    for atom_name, atom in atomgroup.atoms():
    #        num_of_residues = self._get_num_of_residues()
    #        if atom.symbol == 'Na':
    #            res = bridge.AtomGroup(name = 'Na+')
    #            res.set_atom('Na+ ', atom)
    #            chain.set_group(num_of_residues +1, res)
    #        elif atom.symbol == 'Cl':
    #            res = bridge.AtomGroup(name = 'Cl-')
    #            res.set_atom('Cl- ', atom)
    #            chain.set_group(num_of_residues +1, res)
    #        else:
    #            raise
    #    # num_of_chains = self._model.get_number_of_groups()
    #    # self._model.set_group(chr(ord('A') +((num_of_chains +1) % 26)), chain)
    #    self._model.set_group('Z', chain)

    def _add_ions(self, atomgroup, ions):
        assert isinstance(atomgroup, bridge.AtomGroup)
        assert isinstance(ions, bridge.AtomGroup)

        count = 0
        for atom_name, atom in ions.atoms():
            # check collision of the ion index
            new_name = ''
            while True:
                new_name = '{atom_name}{count}'.format(atom_name=atom_name,
                                                       count=count)
                if not atomgroup.has_atomkey(new_name):
                    break
                count += 1
            atomgroup.set_atom(new_name, atom)

    def _get_num_of_residues(self):
        answer = 0
        for chain_name, chain in self._model.groups():
            answer += chain.get_number_of_groups()
        return answer

    def _check_ionpairs_by_dummy(self):
        select_Na = bridge.Select_Atom('Na')
        select_Cl = bridge.Select_Atom('Cl')
        model_Na = self._model.select(select_Na)
        model_Cl = self._model.select(select_Cl)

        # print(model_Na)
        # print(model_Cl)

        atom_list = []
        for chain_name, chain in model_Na.groups():
            for resid, res in chain.groups():
                for atomname, atom in res.atoms():
                    select_range = bridge.Select_Range(atom.xyz, 4.0)
                    ionpairs = model_Cl.select(select_range)
                    if ionpairs.get_number_of_all_atoms() > 0:
                        for chain_name2, chain2 in ionpairs.groups():
                            for resid2, res2 in chain2.groups():
                                for atomname2, atom2 in res2.atoms():
                                    self._logger.info(
                                        'found ion pair: {}'.format(atom.path))
                                    self._logger.info(
                                        '              : {}'.format(atom2.path))
                                    atom_list.append(atom.path)
                                    atom_list.append(atom2.path)

        for atom_path in atom_list:
            (model_name, chain_name, resid, atom_name) = self._divide_path(atom_path)
            self._model[chain_name][resid].erase_atom(atom_name)

    def _check_ionpairs(self):
        self._logger.info("_check_ionpairs()")
        select_Na = bridge.Select_Atom('Na')
        select_Cl = bridge.Select_Atom('Cl')
        model_Na = self._model.select(select_Na)
        model_Cl = self._model.select(select_Cl)

        ip = bridge.IonPair(self._model)
        ionpairs = ip.get_ion_pairs()

        for pair in ionpairs:
            anion_path = pair[0]
            cation_path = pair[1]
            self._logger.info(
                'pair> {} <-> {}'.format(anion_path, cation_path))

    def _reorder_ions_for_amber(self, atomgroup):
        assert isinstance(atomgroup, bridge.AtomGroup)
        self._logger.debug(str(atomgroup))

        def reorder_ions_core(atomgroup):
            ions = []
            for key, subgrp in atomgroup.groups():
                new_ions = reorder_ions_core(subgrp)
                ions.extend(new_ions)
            for key, atom in atomgroup.atoms():
                if atom.symbol not in ('Na', 'Cl'):
                    self._logger.debug(
                        "pass>'{}':'{}'@{}".format(key, atom.symbol, atomgroup.name))
                else:
                    self._logger.debug(
                        "FOUND>'{}':'{}'@{}".format(key, atom.symbol, atomgroup.name))
                    atomgroup.erase_atom(key)
                    ions.append(atom)
            return ions

        new_atomgroup = bridge.AtomGroup(atomgroup)
        ions = reorder_ions_core(new_atomgroup)

        # atomgroup for ions
        chain = brdige.AtomGroup()
        resid = 1
        for ion in ions:
            res = brdige.AtomGroup()
            res.name = ion.name
            res.set_atom(ion.symbol, ion)
            chain.set_group(resid, res)
            resid += 1

        num_of_chains = atomgroup.get_number_of_groups()
        new_chain_id = chr(ord('A') + (num_of_chains % 26))
        chain.name = new_chain_id
        new_atomgroup.set_group(new_chain_id, chain)

        return new_atomgroup

    def _prepare_mdin_step1(self):
        mdin_templ = """
        #
        &cntrl
          imin=1, ntmin=3, maxcyc=50000, ncyc=0, drms=0.0001,
          ntpr=200,
          ntb=0, ntc=1, ntf=1, cut=9999.999,
          igb=1,
        &end
        &lmod
          xmin_method="LBFGS"
        &end
        """
        mdin_templ = mdin_templ.lstrip('\n')
        mdin_templ = brdige.Utils.unindent_block(mdin_templ)
        mdin_templ = brdige.Utils.add_spaces(mdin_templ, 2)

        f = open('md1.in', 'w')
        f.write(mdin_templ)
        f.close()

    def _prepare_mdin_step2(self):
        mdin_templ = """
        #
        &cntrl
          imin=0, irest=0, ntx=1,
          ntt=1, temp0=300.0, tautp=0.1, tempi=0.0,
          ntb=0, ntc=2, ntf=2, cut=10.0,
          ntpr=5000, ntwx=500, ntwe=500,
          nstlim=500000, dt=0.002,
          ibelly=1,
          bellymask='${belly_mask}'
        &end
        """
        mdin_templ = mdin_templ.lstrip('\n')
        mdin_templ = brdige.Utils.unindent_block(mdin_templ)
        mdin_templ = brdige.Utils.add_spaces(mdin_templ, 2)

        (start, end) = self._get_wat_residues()
        belly_mask = ':{start}-{end}'.format(start=start,
                                             end=end)

        mdin_contents = string.Template(mdin_templ).substitute({
            'belly_mask': belly_mask
        })

        f = open('md2.in', 'w')
        f.write(mdin_contents)
        f.close()

    def _prepare_mdin_step3(self):
        mdin_templ = """
        #
        # relax
        &cntrl
          imin=0, ntmin=3, maxcyc=50000, ncyc=0, drms=0.0001,
          ntpr=200,
          ntb=0, ntc=1, ntf=1, cut=9999.999,
          ibelly=1,
          bellymask='${belly_mask}'
        &end
        &lmod
          xmin_method="LBFGS"
        &end
        """
        mdin_templ = mdin_templ.lstrip('\n')
        mdin_templ = brdige.Utils.unindent_block(mdin_templ)
        mdin_templ = brdige.Utils.add_spaces(mdin_templ, 2)

        (start, end) = self._get_ion_residues()
        belly_mask = ':{start}-{end}'.format(start=start,
                                             end=end)

        mdin_contents = string.Template(mdin_templ).substitute({
            'belly_mask': belly_mask
        })

        f = open('md3.in', 'w')
        f.write(mdin_contents)
        f.close()

    def _get_wat_residues(self):
        # load PDB file
        pdb_obj = brdige.Pdb()
        # pdb_obj.debug = debug
        pdb_filepath = 'md2_before.pdb'
        pdb_obj.load(pdb_filepath)
        models = pdb_obj.get_atomgroup()

        start = 999999
        end = -1

        select_WAT = brdige.Select_Name('WAT')
        WATs = models.select(select_WAT)
        for model_name, model in WATs.groups():
            for chain_name, chain in model.groups():
                for resid, res in chain.groups():
                    i = int(resid)
                    start = min(start, i)
                    end = max(end, i)

        return (start, end)

    def _get_ion_residues(self, pdb_filepath='md3_before.pdb'):
        # load PDB file
        pdb_obj = brdige.Pdb()
        # pdb_obj.debug = debug
        pdb_filepath = 'md3_before.pdb'
        pdb_obj.load(pdb_filepath)
        models = pdb_obj.get_atomgroup()

        start = 999999
        end = -1

        select_Na = brdige.Select_Name('Na+')
        Na_s = models.select(select_Na)
        for model_name, model in Na_s.groups():
            for chain_name, chain in model.groups():
                for resid, res in chain.groups():
                    i = int(resid)
                    start = min(start, i)
                    end = max(end, i)

        select_Cl = brdige.Select_Name('Cl-')
        Cl_s = models.select(select_Na)
        for model_name, model in Cl_s.groups():
            for chain_name, chain in model.groups():
                for resid, res in chain.groups():
                    start = min(start, resid)
                    end = max(end, resid)

        return (start, end)

    def _do_leap(self, leapin_path):
        cmdline = 'tleap -s -f {}'.format(leapin_path)
        self._exec_cmd(cmdline)

    def _make_prep(self, ligand):
        # setup ligand_charge_table
        ligand_charge_table = {}
        ligand_charge_table['FAD'] = -2

        self._logger.info('>>>> make prep: {}'.format(ligand))

        # select ligand
        ligand_selecter = brdige.Select_Name(ligand)
        ag_ligand = brdige.AtomGroup()
        ag_ligand[ligand] = self._model.select(ligand_selecter)

        # make ligand pdb
        biopdb = brdige.Pdb(mode='amber')
        biopdb.set_by_atomgroup(ag_ligand)
        ligand_pdb_filepath = '{ligand}.pdb'.format(ligand=ligand)
        self._logger.info(
            'save ligand pdb file: {}'.format(ligand_pdb_filepath))
        pdb_f = open(ligand_pdb_filepath, 'w')
        pdb_f.write(str(biopdb))
        pdb_f.close()

        #
        self._logger.info('calc charges...')
        charge = 0
        if ligand in ligand_charge_table:
            charge = ligand_charge_table[ligand]

        # exec antechamber
        antechamber_cmd = '{ANTECHAMBER} -fi pdb -i {LIG}.pdb -fo prepi -o {LIG}.prep -at gaff -c bcc -nc {CHARGE} -rn {LIG}'.format(
            ANTECHAMBER=self._antechamber_cmd,
            LIG=ligand,
            CHARGE=charge
        )
        self._logger.debug('run: {cmd}'.format(cmd=antechamber_cmd))
        self._exec_cmd(antechamber_cmd)

        #
        self._logger.info('parmchk...')
        parmchk_cmd = '{PARMCHK} -i {LIG}.prep -f prepi -o {LIG}.frcmod'.format(
            PARMCHK=self._parmchk_cmd,
            LIG=ligand
        )
        self._logger.debug('run: {cmd}'.format(cmd=parmchk_cmd))
        self._exec_cmd(parmchk_cmd)

        self._logger.info('<<<< make prep: {} done.'.format(ligand))

    def _exec_cmd(self, cmdline):
        self._logger.debug('exec cmd: "{}"'.format(cmdline))
        cmd = shlex.split(cmdline)
        proc = subprocess.Popen(cmd,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        retcode = proc.returncode

        if retcode != 0:
            output = stdout
            output += '\n'
            output += stderr
            raise qclo.QcScriptRunningError(cmd, output)
        else:
            return True

    def _make_group_DB(self):
        """
        leapで処理できる残基名(?)をregistered_groupに取得
        """
        leapin_templ = """
        ${source_leaprc}
        quit
        """
        leapin_templ = leapin_templ.lstrip('\n')
        leapin_templ = brdige.Utils.unindent_block(leapin_templ)

        leapin_contents = string.Template(leapin_templ).substitute({
            'source_leaprc': self._get_leaprc_source_lines()
        })
        fout = open(self._leapin_premd_filepath, 'w')
        fout.write(leapin_contents)
        fout.close()

        self._do_leap(self._leapin_premd_filepath)

        self._registered_group = set()
        re_grp_block = re.compile('^Loading:\s+(\S+)')
        f_leaplog = open('leap.log', 'r')
        for line in f_leaplog:
            line = line.strip()
            match_obj = re_grp_block.match(line)
            if match_obj != None:
                grp = match_obj.group(1)
                self._registered_group.add(grp)
                self._logger.debug(' registered group: {}'.format(grp))
        f_leaplog.close()

        # add 'WAT'
        self._registered_group.add('WAT')

    def _make_res_DB(self):
        """
        """
        biopdb = brdige.Pdb(mode='amber')
        protein = brdige.AtomGroup()
        protein.set_group('model_1', self._model)
        protein_amber = biopdb.get_modpdb_atomgroup(protein)

        self._res_set = set()
        for chain_key, chain in protein_amber['model_1'].groups():
            for res_key, res in chain.groups():
                self._res_set.add(res.name)

        self._logger.debug('input group:> ')
        for grp in self._res_set:
            self._logger.info(' {}'.format(grp))


if __name__ == "__main__":
    main()
