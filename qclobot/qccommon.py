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

import os
import tempfile

def get_tmpfile_path(suffix='', prefix='tmp', tmp_dir=''):
    '''
    return tmpfile path
    '''
    if tmp_dir == '':
        if 'TEMP' in os.environ:
            tmp_dir = os.environ['TEMP']
        elif 'TMP' in os.environ:
            tmp_dir = os.environ['TMP']
        else:
            tmp_dir = '/tmp'

    try:
        (tmpfile_h, tmpfile_path) = tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=tmp_dir)
    except:
        print(suffix, prefix, tmp_dir)
        raise
    os.close(tmpfile_h)
    return tmpfile_path
    

    
