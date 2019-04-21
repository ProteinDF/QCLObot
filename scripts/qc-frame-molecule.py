#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import logging
import logging.config

import pprint

try:
    import msgpack
except:
    import msgpack_pure as msgpack

import proteindf_bridge as bridge
import proteindf_tools as pdf
import qclobot as qclo


def main():
    parser = argparse.ArgumentParser(description='create QCLO frame PDB file')
    parser.add_argument('frame',
                        nargs=1,
                        help='frame directory')
    parser.add_argument('-o', '--output_path',
                        nargs=1,
                        default=[""])
    parser.add_argument('-c', '--calc_charges_at_iteration',
                        nargs='?',
                        default=0,
                        const=-1,
                        help='calculate and store charges at specified iteration')
    parser.add_argument('-r', '--reference',
                        nargs=1,
                        default=[""])
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    frame_name = args.frame[0]
    reference_path = args.reference[0]
    output_path = args.output_path[0]
    #if len(output_path) == 0:
    #    output_path = "{frame_name}.brd".format(frame_name=frame_name)
    calc_charge = args.calc_charges_at_iteration
    verbose = args.verbose

    if verbose:
        print("frame: {}".format(frame_name))

    # frame
    if os.path.isdir(frame_name):
        frame = qclo.QcFrame(frame_name)
        frame_molecule = bridge.AtomGroup(frame.frame_molecule)

        if calc_charge != 0:
            if verbose:
                print("calculate Mulliken population")

            pop_vtr = frame.pop(iteration = calc_charge)
            frame_molecule.assign_charges(pop_vtr)

        if reference_path:
            reference = bridge.load_atomgroup(reference_path)
            frame_molecule = frame_molecule.restructure(reference)

        if output_path:
            if verbose:
                print("output brd file: {}".format(output_path))
            bridge.save_atomgroup(frame_molecule, output_path)
        else:
            print(frame_molecule)

    else:
        print("{} is not directory.".format(frame_name))


if __name__ == '__main__':
    main()
