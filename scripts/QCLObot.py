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

import argparse
import logging
import logging.config

try:
    import msgpack
except:
    import msgpack_pure as msgpack

import pdfbridge as bridge
import pdfpytools as pdf
import qclobot as qclo

def main():
    parser = argparse.ArgumentParser(description='QCLObot: QM solver based on QCLO method for large-system.')
    parser.add_argument('senario_file_path',
                        nargs=1,
                        help='QCLO senario file (YAML_format)')
    parser.add_argument('-l', '--logconfig',
                        nargs=1,
                        action='store',
                        help='logconfig file')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    senario_file_path = args.senario_file_path[0]
    logconfig_path = ''
    if args.logconfig:
        logconfig_path = args.logconfig[0]

    setup_logconfig(logconfig_path)

    print(senario_file_path)
    qcctrl = qclo.QcControl()
    qcctrl.run(senario_file_path)
    
        
def setup_logconfig(logconfig_path):
    if logconfig_path:
        logging.config.fileConfig(logconfig_path)
    else:
        logging.basicConfig(
            level=logging.WARNING,
            format='%(asctime)s[%(levelname)s]%(name)s %(message)s'
        )
        logfile = logging.FileHandler('qclobot.log')
        logfile.setLevel(logging.INFO)
        logging.getLogger().addHandler(logfile)

        console = logging.StreamHandler()
        console.setLevel(logging.WARNING)
        logging.getLogger().addHandler(console)

        
def load_brdfile(brdfile_path):
    brdfile = open(brdfile_path, 'rb')
    brddata = msgpack.unpackb(brdfile.read())
    brdfile.close()
    atomgroup = bridge.AtomGroup(brddata)

    return atomgroup

    
if __name__ == '__main__':
    main()
    
    #import cProfile
    #cProfile.run('main()', 'qclo_prof')

    #import pstats
    #ps = pstats.Stats('qclo_prof')
    #ps.sort_stats('cumulative').print_stats()
    #ps.sort_stats('time').print_stats()

