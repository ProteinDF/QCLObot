#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import logging
logger = logging.getLogger(__name__)

import pdfbridge
from .taskobject import TaskObject

class QcNeutralize(TaskObject):
    ''' neutralize protein
    '''
    def __init__(self, *args, **kwargs):
        super(QcNeutralize, self).__init__(*args, **kwargs)

    # ==================================================================
    # properties
    # ==================================================================
    def _get_input_pdb_filepath(self):
        path = self._data.get('input_pdb_filepath', 'input.pdb')
        return path
    input_pdb_filepath = property(_get_input_pdb_filepath)


    def _get_output_pdb_filepath(self):
        path = self._data.get('output_pdb_filepath', 'output.pdb')
        return path
    output_pdb_filepath = property(_get_output_pdb_filepath)

    # ==================================================================
    # method
    # ==================================================================
    def run(self):
        neutral_model = self._neutralize(self.model)
        self.output_model = self._reorder_ions_for_amber(neutral_model)

        return self


    def _neutralize(self, model):
        ip = pdfbridge.IonPair(model)
        ionpairs = ip.get_ion_pairs()

        # 処理しやすいように並べ替え
        exempt_list = [] # 免除リスト
        for (anion_path, cation_path, anion_type, cation_type) in ionpairs:
            (anion_chain_name, anion_res_name) = pdfbridge.AtomGroup.divide_path(anion_path)
            (cation_chain_name, cation_res_name) = pdfbridge.AtomGroup.divide_path(cation_path)
            exempt_list.append((anion_chain_name, anion_res_name, anion_type))
            exempt_list.append((cation_chain_name, cation_res_name, cation_type))
            logger.info("ion pair found: {anion_path} <-> {cation_path}".format(
                anion_path=anion_path,
                cation_path=cation_path)
            )

        output_model = pdfbridge.AtomGroup(model)
        modeling = pdfbridge.Modeling()
        for chain_name, chain in output_model.groups():
            for resid, res in chain.groups():
                resname = res.name

                if res.has_atom('H3'):
                    if (chain_name, resid, 'NTM') not in exempt_list:
                        # N-term
                        ag = modeling.neutralize_Nterm(res)
                        logger.info("add ion for N-term: {}".format(ag))
                        self._add_ions(res, ag)
                    else:
                        logger.info('exempt adding ion: {}/{} Nterm'.format(chain_name, resname))
                if res.has_atom('OXT'):
                    if (chain_name, resid, 'CTM') not in exempt_list:
                        # C-term
                        ag = modeling.neutralize_Cterm(res)
                        logger.info("add ion for C-term: {}".format(ag))
                        self._add_ions(res, ag)
                    else:
                        logger.info('exempt adding ion: {}/{} Cterm'.format(chain_name, resname))

                if resname == 'GLU':
                    if (chain_name, resid, 'GLU') not in exempt_list:
                        ag = modeling.neutralize_GLU(res)
                        logger.info("add ion for GLU({}): {}".format(resid, ag))
                        self._add_ions(res, ag)
                    else:
                        logger.info('exempt adding ion: {}/{} GLU'.format(chain_name, resname))
                elif resname == 'ASP':
                    if (chain_name, resid, 'ASP') not in exempt_list:
                        ag = modeling.neutralize_ASP(res)
                        logger.info("add ion for ASP({}): {}".format(resid, ag))
                        self._add_ions(res, ag)
                    else:
                        logger.info('exempt adding ion: {}/{} ASP'.format(chain_name, resname))
                elif resname == 'LYS':
                    if (chain_name, resid, 'LYS') not in exempt_list:
                        ag = modeling.neutralize_LYS(res)
                        logger.info("add ion for LYS({}): {}".format(resid, ag))
                        self._add_ions(res, ag)
                    else:
                        logger.info('exempt adding ion: {}/{} LYS'.format(chain_name, resname))
                elif resname == 'ARG':
                    if (((chain_name, resid, 'ARG') not in exempt_list) and
                        ((chain_name, resid, 'ARG1') not in exempt_list) and
                        ((chain_name, resid, 'ARG2') not in exempt_list)):
                        ag = modeling.neutralize_ARG(res)
                        logger.info("add ion for ARG({}): {}".format(resid, ag))
                        self._add_ions(res, ag)
                    else:
                        logger.info('exempt adding ion: {}/{} ARG'.format(chain_name, resname))

                elif resname == 'FAD':
                    ag = modeling.neutralize_FAD(res)
                    logger.info("add ion for FAD({}): {}".format(resid, ag))
                    self._add_ions(res, ag)

        return output_model


    def _add_ions(self, atomgroup, ions):
        assert isinstance(atomgroup, pdfbridge.AtomGroup)
        assert isinstance(ions, pdfbridge.AtomGroup)

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


    def _reorder_ions_for_amber(self, model):
        assert isinstance(model, pdfbridge.AtomGroup)

        def pick_ions(atomgroup):
            ions = []
            for key, subgrp in atomgroup.groups():
                new_ions = pick_ions(subgrp)
                ions.extend(new_ions)
            for key, atom in atomgroup.atoms():
                if atom.symbol not in ('Na', 'Cl'):
                    logger.debug("pass>'{}':'{}'@{}".format(key, atom.symbol, atomgroup.name))
                else:
                    logger.debug("FOUND>'{}':'{}'@{}".format(key, atom.symbol, atomgroup.name))
                    atomgroup.erase_atom(key)
                    ions.append(atom)
            return ions

        new_model = pdfbridge.AtomGroup(model)
        ions = pick_ions(new_model)

        # atomgroup for ions
        chain = pdfbridge.AtomGroup()
        resid = 1
        for ion in ions:
            res = pdfbridge.AtomGroup()
            res.name = ion.name
            res.set_atom(ion.symbol, ion)
            chain.set_group(resid, res)
            resid += 1

        new_chain_id = chr(ord(self._get_last_chainid(new_model)) + 1)
        chain.name = new_chain_id
        new_model.set_group(new_chain_id, chain)

        return new_model

    def _get_last_chainid(self, model):
        max_chain_id = 0
        for chain_id, chain in model.groups():
            max_chain_id = max(ord(chain_id) - ord('A'), max_chain_id)
        last_chainid = chr(ord('A') + max_chain_id)
        return last_chainid
