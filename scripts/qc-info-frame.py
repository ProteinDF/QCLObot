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
    parser = argparse.ArgumentParser(description='QCLO frame information')
    parser.add_argument('frame',
                        nargs=1,
                        help='frame directory')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    # print(args)

    frame_name = args.frame[0]
    print("frame: {}".format(frame_name))

    if os.path.isdir(frame_name):
        frame = qclo.QcFrame(frame_name)
        # print(frame)

        num_of_AOs = frame.get_number_of_AOs()
        print("AOs: {}".format(num_of_AOs))

        frame_molecule = frame.frame_molecule
        # print(frame_molecule)

        is_finished_scf = frame.is_finished_scf
        print("SCF: {}".format(is_finished_scf))
        if is_finished_scf:
            #summary = frame.summary()
            #print(summary)
            pop_vtr = frame.pop()
            # print(pop_vtr)

            assign_charges(frame_molecule, pop_vtr)
            print(frame_molecule)


    else:
        print("{} is not directory.".format(frame))


if __name__ == '__main__':
    import cProfile

    pr = cProfile.Profile()
    pr.enable()

    main()

    pr.disable()
    pr.dump_stats('info-frame.profile')
