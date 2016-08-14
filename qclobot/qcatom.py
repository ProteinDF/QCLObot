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

import pdfbridge as bridge
import pdfpytools as pdf

logger = logging.getLogger(__name__)


class QcAtom(bridge.Atom):
    def __init__(self, *args, **kwargs):
        super(QcAtom, self).__init__(*args, **kwargs)

        if kwargs.get('debug'):
            logger.addHandler(logging.StreamHandler())
            logger.setLevel(logging.DEBUG)
        else:
            logger.addHandler(logging.NullHandler())
            logger.setLevel(logging.INFO)

        suffix = "." + self.symbol
        self.basisset = kwargs.get('basisset', 'O-DZVP2' + suffix)
        self.basisset_j = kwargs.get('basisset_j', 'A-DZVP2' + suffix)
        self.basisset_xc = kwargs.get('basisset_xc', 'A-DZVP2' + suffix)
        self.basisset_gridfree = kwargs.get('basisset_gridfree', 'O-DZVP2' + suffix)

        if isinstance(args[0], QcAtom):
            self.basisset = args[0].basisset
            self.basisset_j = args[0].basisset_j
            self.basisset_xc = args[0].basisset_xc
            self.basisset_gridfree = args[0].basisset_gridfree
        
        self._qc_parent = kwargs.get('qc_parent', None)

    # transform to bridge.Atom
    def get_Atom(self):
        atm = bridge.Atom(self)
        return atm

    def set_by_raw_data(self, data):
        super(QcAtom, self).set_by_raw_data(data)

        suffix = "." + self.symbol
        self.basisset = data.get('basisset', 'O-DZVP2' + suffix)
        self.basisset_j = data.get('basisset_j', 'A-DZVP2' + suffix)
        self.basisset_xc = data.get('basisset_xc', 'A-DZVP2' + suffix)
        self.basisset_gridfree = data.get('basisset_gridfree', 'O-DZVP2' + suffix)
        
    def get_raw_data(self):
        data = super(QcAtom, self).get_raw_data()
        data['basisset'] = self.basisset
        data['basisset_j'] = self.basisset_j
        data['basisset_xc'] = self.basisset_xc
        data['basisset_gridfree'] = self.basisset_gridfree

        return data
        
    #
    def get_number_of_AOs(self):
        basis2 = pdf.Basis2()
        bs = basis2.get_basisset(self.basisset)
        return bs.get_number_of_AOs()
        
    # basisset -------------------------------------------------------
    def _get_basisset(self):
        return self._basisset

    def _set_basisset(self, name):
        self._basisset = bridge.Utils.to_unicode(name)

    basisset = property(_get_basisset, _set_basisset)

    # basisset_j ------------------------------------------------------
    def _get_basisset_j(self):
        return self._basisset_j

    def _set_basisset_j(self, name):
        self._basisset_j = bridge.Utils.to_unicode(name)

    basisset_j = property(_get_basisset_j, _set_basisset_j)

    # basisset_xc -----------------------------------------------------
    def _get_basisset_xc(self):
        return self._basisset_xc

    def _set_basisset_xc(self, name):
        self._basisset_xc = bridge.Utils.to_unicode(name)

    basisset_xc = property(_get_basisset_xc, _set_basisset_xc)

    # basisset_gridfree -----------------------------------------------
    def _get_basisset_gridfree(self):
        return self._basisset_gridfree

    def _set_basisset_gridfree(self, name):
        self._basisset_gridfree = bridge.Utils.to_unicode(name)

    basisset_gridfree = property(_get_basisset_gridfree, _set_basisset_gridfree)

    # atom label ------------------------------------------------------
    def _get_atomlabel(self):
        symbol = self.symbol
        label = self.label
        kind = symbol
        if len(label) > 0:
            kind = '{}@{}'.format(symbol, label)
        return kind

    atomlabel = property(_get_atomlabel)



