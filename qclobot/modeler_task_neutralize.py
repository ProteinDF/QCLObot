#!/usr/bin/env python
# -*- coding: utf-8 -*-

import proteindf_bridge as bridge

from .modeler_task_object import ModelerTaskObject

import logging
logger = logging.getLogger(__name__)


class ModelerTaskNeutralize(ModelerTaskObject):
    ''' neutralize protein
    '''
    _task_name = "neutralize"

    def __init__(self, parent, task):
        super().__init__(parent, task)

    # ==================================================================
    # properties
    # ==================================================================

    # ==================================================================
    # method
    # ==================================================================
    def run(self):
        neutral_model = self._neutralize(self.model)
        self.output_model = self._reorder_ions_for_amber(neutral_model)

        return self

    def _neutralize(self, model):
        ip = bridge.IonPair(model)
        ionpairs = ip.get_ion_pairs()

        # 処理しやすいように並べ替え
        exempt_list = []  # 免除リスト
        for (anion_path, cation_path, anion_type, cation_type) in ionpairs:
            (anion_chain_name, anion_res_name) = bridge.AtomGroup.divide_path(anion_path)
            (cation_chain_name, cation_res_name) = bridge.AtomGroup.divide_path(cation_path)
            exempt_list.append((anion_chain_name, anion_res_name, anion_type))
            exempt_list.append((cation_chain_name, cation_res_name, cation_type))
            logger.info("ion pair found: {anion_path} <-> {cation_path}".format(
                anion_path=anion_path,
                cation_path=cation_path)
            )

        output_model = bridge.AtomGroup(model)
        modeling = bridge.Modeling()
        for chain_name, chain in output_model.groups():
            for resid, res in chain.groups():
                resname = res.name

                if resname in ['ACE', 'NME']:
                    continue

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
                    if (((chain_name, resid, 'ARG') not in exempt_list) and ((chain_name, resid, 'ARG1') not in exempt_list) and ((chain_name, resid, 'ARG2') not in exempt_list)):
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

    def _reorder_ions_for_amber(self, model):
        logger.info("reorder ions: start")
        assert isinstance(model, bridge.AtomGroup)

        def pick_ions(atomgroup):
            ions = []
            for key, subgrp in atomgroup.groups():
                new_ions = pick_ions(subgrp)
                ions.extend(new_ions)

            remove_atom_keys = []
            for key, atom in atomgroup.atoms():
                if atom.symbol not in ('Na', 'Cl'):
                    logger.debug("pass>'{}':'{}'@{}".format(key, atom.symbol, atomgroup.name))
                else:
                    logger.debug("FOUND>'{}':'{}'@{}".format(key, atom.symbol, atomgroup.name))
                    ions.append(atom)
                    remove_atom_keys.append(key)

            for key in remove_atom_keys:
                atomgroup.remove_atom(key)

            return ions

        new_model = bridge.AtomGroup(model)
        ions = pick_ions(new_model)
        logger.info("ions: {}".format(len(ions)))

        # atomgroup for ions
        chain = bridge.AtomGroup()
        resid = 1
        for ion in ions:
            res = bridge.AtomGroup()
            res.name = ion.name
            res.set_atom(ion.symbol, ion)
            chain.set_group(resid, res)
            resid += 1

        new_chain_id = chr(ord(self._get_last_chainid(new_model)) + 1)
        chain.name = new_chain_id
        new_model.set_group(new_chain_id, chain)

        logger.info("reorder ions: end")
        return new_model

    def _get_last_chainid(self, model):
        max_chain_id = 0
        for chain_id, chain in model.groups():
            max_chain_id = max(ord(chain_id) - ord('A'), max_chain_id)
        last_chainid = chr(ord('A') + max_chain_id)
        return last_chainid
