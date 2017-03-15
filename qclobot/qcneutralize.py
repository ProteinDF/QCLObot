#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import logging
logger = logging.getLogger(__name__)

from pdfbridge
from .taskobject import TaskObject

class QcNeutralize(TaskObject):
    ''' neutralize protein
    '''
    def __init__(self, name):
        super(QcNeutralize, self).__init__(name=name)


    def run(self):
        ip = pdfbridge.IonPair(self._model)
        ionpairs = ip.get_ion_pairs()
       
        # 処理しやすいように並べ替え
        exempt_list = []
        for (anion_path, cation_path, anion_type, cation_type) in ionpairs:
            (anion_chain_name, anion_res_name) = pdfbridge.AtomGroup.divide_path(anion_path)
            (cation_chain_name, cation_res_name) = pdfbridge.AtomGroup.divide_path(cation_path)
            exempt_list.append((anion_chain_name, anion_res_name, anion_type))
            exempt_list.append((cation_chain_name, cation_res_name, cation_type))

        output_model = pdfbridge.AtomGroup(self.model)
        modeling = pdfbridge.Modeling()
        for chain_name, chain in output_model.groups():
            for resid, res in chain.groups():
                resname = res.name

                if res.has_atom('H3'):
                    if (chain_name, resid, 'NTM') not in exempt_list:
                        # N-term
                        ag = modeling.neutralize_Nterm(res)
                        self._logger.info("add ion for N-term: {}".format(ag))
                        self._add_ions(res, ag)
                    else:
                        self._logger.info('exempt adding ion: {}/{} Nterm'.format(chain_name, resname))
                if res.has_atom('OXT'):
                    if (chain_name, resid, 'CTM') not in exempt_list:
                        # C-term
                        ag = modeling.neutralize_Cterm(res)
                        self._logger.info("add ion for C-term: {}".format(ag))
                        self._add_ions(res, ag)
                    else:
                        self._logger.info('exempt adding ion: {}/{} Cterm'.format(chain_name, resname))

                if resname == 'GLU':
                    if (chain_name, resid, 'GLU') not in exempt_list:
                        ag = modeling.neutralize_GLU(res)
                        self._logger.info("add ion for GLU({}): {}".format(resid, ag))
                        self._add_ions(res, ag)                        
                    else:
                        self._logger.info('exempt adding ion: {}/{} GLU'.format(chain_name, resname))
                elif resname == 'ASP':
                    if (chain_name, resid, 'ASP') not in exempt_list:
                        ag = modeling.neutralize_ASP(res)
                        self._logger.info("add ion for ASP({}): {}".format(resid, ag))
                        self._add_ions(res, ag)
                    else:
                        self._logger.info('exempt adding ion: {}/{} ASP'.format(chain_name, resname))
                elif resname == 'LYS':
                    if (chain_name, resid, 'LYS') not in exempt_list:
                        ag = modeling.neutralize_LYS(res)
                        self._logger.info("add ion for LYS({}): {}".format(resid, ag))
                        self._add_ions(res, ag)
                    else:
                        self._logger.info('exempt adding ion: {}/{} LYS'.format(chain_name, resname))
                elif resname == 'ARG':
                    if (((chain_name, resid, 'ARG') not in exempt_list) and
                        ((chain_name, resid, 'ARG1') not in exempt_list) and
                        ((chain_name, resid, 'ARG2') not in exempt_list)):
                        ag = modeling.neutralize_ARG(res)
                        self._logger.info("add ion for ARG({}): {}".format(resid, ag))
                        self._add_ions(res, ag)
                    else:
                        self._logger.info('exempt adding ion: {}/{} ARG'.format(chain_name, resname))

                elif resname == 'FAD':
                    ag = modeling.neutralize_FAD(res)
                    self._logger.info("add ion for FAD({}): {}".format(resid, ag))
                    self._add_ions(res, ag)
 
        self.output_model = model


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

        
        
