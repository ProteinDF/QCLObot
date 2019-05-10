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

def get_atomgroup(brd_file_path):
    mpac_file = open(brd_file_path, "rb")
    mpac_data = msgpack.unpackb(mpac_file.read())
    mpac_file.close()
    atomgroup = bridge.AtomGroup(mpac_data)

    return atomgroup

def save_atomgroup(atomgroup, file_path):
    data = atomgroup.get_raw_data()
    mpac = msgpack.packb(data)

    with open(file_path, "wb") as f:
        f.write(mpac)



def main():
    parser = argparse.ArgumentParser(description='create QCLO frame PDB file')
    parser.add_argument('frame',
                        nargs=1,
                        help='frame directory')
    parser.add_argument('ref_brd_path',
                        nargs=1,
                        help='reference brd file')
    parser.add_argument('-o', '--output_path',
                        nargs=1,
                        default=["output.brd"])
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    # print(args)

    frame_name = args.frame[0]
    ref_brd_path = args.ref_brd_path[0]
    output_path = args.output_path[0]
    verbose = args.verbose

    if verbose:
        print("frame: {}".format(frame_name))
        print("reference: {}".format(ref_brd_path))

    # load reference
    ref_ag = get_atomgroup(ref_brd_path)

    # frame
    if os.path.isdir(frame_name):
        frame = qclo.QcFrame(frame_name)

        ag = bridge.AtomGroup(frame.frame_molecule)
        ag_selector = bridge.Select_AtomGroup(ag)

        answer = ref_ag.select(ag_selector)
        # print(answer)

        rest = ref_ag ^ answer
        print(">>>>")
        print(rest)

        if output_path:
            if verbose:
                print("output brd file: {}".format(output_path))
            save_atomgroup(answer, output_path)

    else:
        print("{} is not directory.".format(frame_name))


if __name__ == '__main__':
    import cProfile

    #pr = cProfile.Profile()
    #pr.enable()

    main()

    #pr.disable()
    #pr.dump_stats('program.profile')
