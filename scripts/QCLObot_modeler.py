#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os.path
import logging
import logging.config

import argparse

try:
    import msgpack
except:
    import msgpack_pure as msgpack


import proteindf_bridge as bridge
import proteindf_tools as pdf
import qclobot as qclo

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description='QCLObot: QM solver based on QCLO method for large-system.')
    parser.add_argument('-l', '--logfile',
                        nargs=1,
                        action='store',
                        help='logconfig file')
    parser.add_argument('-L', '--logconfig',
                        nargs=1,
                        action='store',
                        help='logconfig file')
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False)
    parser.add_argument('senario_file_path',
                        nargs=1,
                        help='QCLO senario file (YAML_format)')
    args = parser.parse_args()

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

    logger.info('loading senario: {}'.format(senario_file_path))

    modeler = qclo.QcModeler()
    modeler.run(senario_file_path)


def setup_logging(logfile_path = '', is_debug = False):
    if len(logfile_path) == 0:
        logfile_path = 'qclobot_modeler.log'

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


if __name__ == '__main__':
    if os.path.exists("config.ini"):
        logging.config.fileConfig("config.ini",
                                  disable_exisiting_logger=False)
    main()
