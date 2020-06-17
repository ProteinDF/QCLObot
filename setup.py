#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2019 The ProteinDF development team.
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

import sys
import os
import io
import re

from setuptools import setup
from imp import reload
exec(open("qclobot/_version.py").read())

setup(name='qclobot',
      version=__version__,
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
          'scripts/remove_wat.py',
          'scripts/qc-info-frame.py',
          'scripts/qc-frame-molecule.py',
          'scripts/qc-frame-checkconv.py',
          'scripts/qc-frame-orb.py'
      ],

      install_requires=[
          'configparser',
          'jinja2',
          'rainbow_logging_handler',
          # 'proteindf_bridge',
          # 'proteindf_tools'
      ],
      )
