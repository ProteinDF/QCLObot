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

import pdfbridge as bridge
import pdfpytools as pdf
import qclobot as qclo

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description='QCLObot: QM solver based on QCLO method for large-system.')
    parser.add_argument('senario_file_path',
                        nargs=1,
                        help='QCLO senario file (YAML_format)')
    args = parser.parse_args()
    senario_file_path = args.senario_file_path[0]

    modeler = qclo.QcModeler()
    modeler.run(senario_file_path)


if __name__ == '__main__':
    if os.path.exists("config.ini"):
        logging.config.fileConfig("config.ini",
                                  disable_exisiting_logger=False)
    main()
