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


def assign_charges(atomgroup, charges, charge_index = 0):
    assert(isinstance(atomgroup, bridge.AtomGroup))
    assert(isinstance(charges, bridge.Vector))

    for key, subgrp in atomgroup.groups():
        charge_index = assign_charges(subgrp, charges, charge_index)
    for key, atom in atomgroup.atoms():
        atom.charge = charges.get(charge_index)
        charge_index += 1

    return charge_index


def main():
    parser = argparse.ArgumentParser(description='check convergence in the QCLO frame')
    parser.add_argument('frame',
                        nargs=1,
                        help='frame directory')
    parser.add_argument('-o', '--output_path',
                        nargs=1,
                        default=[""])
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    # print(args)

    frame_name = args.frame[0]
    output_path = args.output_path[0]
    #if len(output_path) == 0:
    #    output_path = "{frame_name}.brd".format(frame_name=frame_name)
    verbose = args.verbose

    if verbose:
        print("frame: {}".format(frame_name))

    # frame
    if os.path.isdir(frame_name):
        frame = qclo.QcFrame(frame_name)
        frame_molecule = bridge.AtomGroup(frame.frame_molecule)

        pop_vtr1 = frame.pop(iteration = 0)
        pop_vtr2 = frame.pop(iteration = -1)
        #print(pop_vtr1)
        #print(pop_vtr2)

        pop_vtr2 -= pop_vtr1
        #print(pop_vtr2)

        frame_molecule.assign_charges(pop_vtr2)

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
