#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import logging
import logging.config

import proteindf_bridge as bridge
import proteindf_tools as pdf
import qclobot as qclo


def main():
    parser = argparse.ArgumentParser(
        description='show frame orbital information')
    parser.add_argument('frame',
                        nargs=1,
                        help='frame directory')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    frame_name = args.frame[0]
    verbose = args.verbose

    if verbose:
        print("frame: {}".format(frame_name))

    if os.path.isdir(frame_name):
        frame = qclo.QcFrame(frame_name)
        orbital_info = frame.get_orbital_info()
        # frame_molecule = bridge.AtomGroup(frame.frame_molecule)
        print(type(orbital_info), len(orbital_info))

        for i, orb in enumerate(orbital_info):
            print("{index:>5}: {orb_str}".format(index=i+1,
                                                 orb_str=str(orb)))


if __name__ == '__main__':
    main()
