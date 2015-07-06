#!/usr/bin/env python

# Copyright (C) 2014 The ProteinDF development team.
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

import sys, os
#from setuptools import setup, find_packages
from distutils.core import setup
from imp import reload

sys.path.append('./qclobot')

setup(name='QCLObot',
      version='2014.0',
      description='building initial guess scripts based on QCLO for the ProteinDF',
      author='Toshiyuki HIRANO',
      author_email='hiracchi@gmail.com',
      url='http://proteindf.github.io/',

      #install_requires = [
      #    'configparser',
      #    'msgpack-python',
      #    'pyyaml',
      #],
      
      packages=['qclobot'],
      scripts=[
          'scripts/qclo_sample.py',
          'scripts/QCLObot.py'
      ],
)
