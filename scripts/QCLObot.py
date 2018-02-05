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

import proteindf_bridge as bridge
import proteindf_tools as pdf
import qclobot as qclo

def main():
    parser = argparse.ArgumentParser(description='QCLObot: QM solver based on QCLO method for large-system.')
    parser.add_argument('senario_file_path',
                        nargs=1,
                        help='QCLO senario file (YAML_format)')
    parser.add_argument('-l', '--logfile',
                        nargs=1,
                        action='store',
                        help='logconfig file')
    parser.add_argument('-L', '--logconfig',
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
    is_debug = args.debug

    senario_file_path = args.senario_file_path[0]

    if args.logconfig:
        # logging module setup by config-file
        logconfig_path = args.logconfig[0]
        logging.config.fileConfig(logconfig_path)
    else:
        # logging module setup
        logfile_path = ''
        if args.logfile:
            logfile_path = args.logfile[0]
        setup_logging(logfile_path, is_debug)

    app_logger = logging.getLogger(__name__)
    app_logger.info('loading senario: {}'.format(senario_file_path))

    qcctrl = qclo.QcControl()
    qcctrl.run(senario_file_path)

    app_logger.info('QCLObot done.')


def setup_logging(logfile_path = '', is_debug = False):
    if len(logfile_path) == 0:
        logfile_path = 'qclobot.log'

    logging_level = logging.INFO
    format_str = '%(asctime)s [%(levelname)s] %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'
    if is_debug:
        logging_level = logging.DEBUG
        format_str ='%(asctime)s [%(levelname)s] [%(name)s] %(message)s'

    logging.basicConfig(
        filename = logfile_path,
        level=logging_level,
        format=format_str,
        datefmt=date_format
    )

    formatter = logging.Formatter(format_str, date_format)

    console = logging.StreamHandler()
    console.setLevel(logging.WARNING)
    console.setFormatter(formatter)
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
