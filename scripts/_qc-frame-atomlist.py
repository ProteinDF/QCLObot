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

def get_fragment_atoms(fragment):
    atoms = []
    for name, subfragment in fragment.groups():
        atoms += get_fragment_atoms(subfragment)

    for name, atom in fragment.atoms():
        atoms += [atom]

    return atoms


def main():
    parser = argparse.ArgumentParser(description='show QCLO frame fragments')
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

    frame_name = args.frame[0]
    output_path = args.output_path[0]
    verbose = args.verbose

    if verbose:
        print("frame: {}".format(frame_name))

    atomgroup = bridge.AtomGroup()
    if not os.path.isdir(frame_name):
        sys.stderr.write("not found: {0}".format(frame_name))

    frame = qclo.QcFrame(frame_name)
    group_index = 0
    for name, fragment in frame.fragments():
        atoms = get_fragment_atoms(fragment)

        subgrp = bridge.AtomGroup()
        for i, atom in enumerate(atoms):
            subgrp.set_atom(i, atom)
        atomgroup.set_group(group_index, subgrp)
        group_index += 1

    if len(output_path) > 0:
        data = atomgroup.get_raw_data()
        if verbose:
            print("output bridge file: {0}".format(output_path))
        bridge.save_msgpack(data, output_path)
    else:
        print(atomgroup)

if __name__ == '__main__':
    main()
