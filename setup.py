#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
import io
import re

from setuptools import setup
from imp import reload
from qclobot import __version__

sys.path.append('./qclobot')

with io.open('qclobot/__init__.py', 'rt', encoding='utf8') as f:
    version = re.search(r'__version__ = \'(.*?)\'', f.read()).group(1)

setup(name='qclobot',
      version=version,
      description='building initial guess scripts based on QCLO for the ProteinDF',
      author='Toshiyuki HIRANO',
      author_email='hiracchi@gmail.com',
      url='http://proteindf.github.io/',
      license='GPLv3',
      packages=['qclobot'],
      zip_safe=False,
      scripts=[
          'scripts/qclo_sample.py',
          'scripts/QCLObot.py',
          'scripts/QCLObot_opt.py',
          'scripts/QCLObot_modeler.py',
          'scripts/relax_protein.py',
          'scripts/relax_protein.sh',
          'scripts/remove_wat.py'
      ],

      install_requires = [
          'configparser',
          'msgpack-python',
          'pyyaml',
          'jinja2',
          #'proteindf_bridge',
          #'proteindf_tools'
      ],
)
