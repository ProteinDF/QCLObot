#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import pprint

import proteindf_bridge as bridge
import proteindf_tools as pdf
import qclobot as qclo

import logging
import logging.config


def main():
    parser = argparse.ArgumentParser(description='apply amber charges')
    parser.add_argument('FILE',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('-o', '--output',
                        nargs='?',
                        help='output bridge file')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    # print(args)

    # setting
    verbose = args.verbose
    mpac_file_path = args.FILE[0]
    output_path = args.output

    # reading
    if (verbose == True):
        print("reading: %s\n" % (mpac_file_path))
    mpac_data = bridge.load_msgpack(mpac_file_path)

    # prepare atomgroup
    atom_group = bridge.AtomGroup(mpac_data)
    model = atom_group.get_group('model_1')
    # print(atom_group)

    # ssbond = bridge.SSBond()
    # ssbond.check(model)
    # print(ssbond.ssbonds)

    # amber
    amber = qclo.AmberObject(name="QM-charge")
    amber.model = model
    amber.charge()

    # output file
    protein = bridge.AtomGroup()
    protein.set_group("model_1", amber.output_model)

    # raw_data = amber.output_model.get_raw_data()
    raw_data = protein.get_raw_data()
    if (verbose == True):
        print("writing: %s\n" % (output_path))
    bridge.save_msgpack(raw_data, output_path)


if __name__ == '__main__':
    main()
