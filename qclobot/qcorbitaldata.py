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

import logging

class QcOrbitalData(object):
    def __init__(self, atom, basisset_name, CGTO_index, basis_type, *args, **kwargs):
        self._logger = logging.getLogger(__name__)
        if kwargs.get('debug'):
            self._logger.addHandler(logging.StreamHandler())
            self._logger.setLevel(logging.DEBUG)
        else:
            self._logger.addHandler(logging.NullHandler())
            self._logger.setLevel(logging.INFO)

        self._atom = atom
        self._basisset_name = basisset_name
        self._CGTO_index = CGTO_index
        self._basis_type = basis_type

    # atom -------------------------------------------------------------
    def _get_atom(self):
        return self._atom
    
    atom = property(_get_atom)
    
    # basisset_name ----------------------------------------------------
    def _get_basisset_name(self):
        return self._basisset_name
        
    basisset_name = property(_get_basisset_name)

    # CGTO_index -------------------------------------------------------
    def _get_CGTO_index(self):
        return self._CGTO_index

    CGTO_index = property(_get_CGTO_index)

    # basis_type -------------------------------------------------------
    def _get_basis_type(self):
        return self._basis_type

    basis_type = property(_get_basis_type)
    
    # operator== -------------------------------------------------------
    def __eq__(self, rhs):
        assert(isinstance(rhs, QcOrbitalData))
        return ((self.atom == rhs.atom) and
                (self.basisset_name == rhs.basisset_name) and
                (self.CGTO_index == rhs.CGTO_index) and
                (self.basis_type == rhs.basis_type))

    # operator!= -------------------------------------------------------
    def __ne__(self, rhs):
        return not self.__eq__(rhs)
    
    
class QcOrbitalInfo(object):
    def __init__(self):
        pass

    def get_orb_info(self):
        basis2 = pdf.Basis2()
        orb_info = []
        
        for subgrp_key, subgrp in qc_fragment.groups():
            orb_info.extend(subgrp.get_orb_info())
        
        for atom_key, qc_atom in qc_fragment.atoms():
            atom = QcAtom(qc_atom)
            bsname = qc_atom.basisset
            basisset = basis2.get_basisset(bsname)
            num_of_CGTOs = len(basisset)
            for CGTO_index in range(num_of_CGTOs):
                CGTO = basisset[CGTO_index]
                shell_type = CGTO.shell_type
                shell_type_id = pdf.ContractedGTO.get_shell_type_id(shell_type)
                num_of_basis_type = shell_type_id * 2 + 1
                for basis_type in range(num_of_basis_type):
                    data = QcOrbitalData(atom = atom,
                                         basisset_name = basisset_name,
                                         CGTO_index = CGTO_index,
                                         basis_type = basis_type)
                    orb_info.append(data)

        return orb_info
            
