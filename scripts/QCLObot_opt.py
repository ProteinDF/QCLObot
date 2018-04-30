#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
    parser = argparse.ArgumentParser(description='QCLObot_opt: QM solver based on QCLO method for large-system.')
    parser.add_argument('bridge_path',
                        nargs=1,
                        help='molecule file (ProteinDF bridge format)')
    parser.add_argument('-t', '--template',
                        nargs=1,
                        help='QCLO senario template file (YAML_format)')
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

    brd_path = args.bridge_path[0]

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

    template_path = None
    if args.template:
        template_path = args.template[0]

    app_logger = logging.getLogger(__name__)

    print('brd_path: {}'.format(brd_path))
    print('template: {}'.format(template_path))
    qcopt = qclo.QcOpt(brd_path = brd_path,
                       template_path = template_path)
    qcopt.run()

    app_logger.info('done.')


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


def get_atom_group(brd_file_path):
    brdfile = open(brd_file_path, 'rb')
    brddata = msgpack.unpackb(brdfile.read())
    brdfile.close()
    atom_group = bridge.AtomGroup(brddata)

    return atom_group


if __name__ == '__main__':
    main()
