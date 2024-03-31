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

import sys
import argparse
import yaml
import rich.logging
import rich.progress


import proteindf_bridge as bridge
import proteindf_tools as pdf
import qclobot as qclo

import logging
import logging.config
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="QCLObot: QM solver based on QCLO method for large-system."
    )
    parser.add_argument(
        "scenario_file_path", nargs=1, help="QCLO scenario file (YAML_format)"
    )
    parser.add_argument(
        "-o", "--output", nargs=1, action="store", help="output log file"
    )
    parser.add_argument(
        "--logconfig", nargs=1, action="store", help="logconfig YAML file"
    )
    # parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("--debug", action="store_true", default=False)
    args = parser.parse_args()

    # setting
    # verbose = args.verbose
    is_debug = args.debug

    scenario_file_path = args.scenario_file_path[0]

    # setup logger
    config_path = ""
    if args.logconfig:
        config_path = args.logconfig[0]

    logfile_path = "qclobot.log"
    if args.output:
        logfile_path = args.output[0]

    setup_logger(is_debug, logfile_path, config_path=config_path)


    # main
    logger.info("loading scenario: {}".format(scenario_file_path))

    qcctrl = qclo.QcControl()
    qcctrl.run(scenario_file_path)

    logger.info("QCLObot done.")


def setup_logger(is_debug=False, logfile_path="", config_path=""):
    if len(config_path) > 0:
        # setup logger by config file
        with open(config_path, 'r') as yml:
            logging_dict = yaml.load(yml, Loader=yaml.SafeLoader)
            logging.config.dictConfig(logging_dict)
    else:
        # setup default logger 
        handlers = []
        level = logging.INFO
        if is_debug:
            level = logging.DEBUG

        rich_handler = rich.logging.RichHandler(rich_tracebacks=True)
        rich_handler.setLevel(level)
        rich_handler.setFormatter(logging.Formatter(fmt="%(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
        handlers.append(rich_handler)

        if len(logfile_path) > 0:
            file_handler = logging.FileHandler(logfile_path)
            file_handler.setLevel(level)
            file_handler.setFormatter(logging.Formatter(fmt="%(asctime)s@%(name)s[%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
            handlers.append(file_handler)

        logging.basicConfig(level=logging.NOTSET, handlers=handlers)


def load_brdfile(brdfile_path):
    atomgroup = bridge.load_atomgroup(brdfile_path)

    return atomgroup


if __name__ == "__main__":
    main()

    # import cProfile
    # cProfile.run('main()', 'qclo_prof')

    # import pstats
    # ps = pstats.Stats('qclo_prof')
    # ps.sort_stats('cumulative').print_stats()
    # ps.sort_stats('time').print_stats()
