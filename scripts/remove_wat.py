#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import msgpack

import proteindf_bridge as bridge


def main():
    # parse args
    parser = argparse.ArgumentParser(description='remove water molecules in bridge file')
    parser.add_argument('FILE',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('-o', '--output',
                        nargs='?',
                        help='output bridge file')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default = False)
    args = parser.parse_args()

    # setting
    mpac_file_path = args.FILE[0]
    output = args.output
    verbose = args.verbose

    # reading
    if verbose:
        print("reading: %s\n" % (mpac_file_path))
    mpac_file = open(mpac_file_path, "rb")
    mpac_data = msgpack.unpackb(mpac_file.read())
    mpac_file.close()

    # prepare atomgroup
    answer = brdige.AtomGroup(mpac_data)
    #print(atomgroup)

    # search 'HOH'
    HOH_selecter = brdige.Select_Name('HOH')
    HOH_grp = answer.select(HOH_selecter)
    answer ^= HOH_grp

    # search 'WAT'
    WAT_selecter = brdige.Select_Name('WAT')
    WAT_grp = answer.select(WAT_selecter)
    answer ^= WAT_grp

    # output
    if output:
        if (verbose == True):
            print("writing: %s\n" % (output))
        savefile = open(output, "wb");
        savefile.write(msgpack.packb(answer.get_raw_data()))
        savefile.close()
    else:
        print(answer)

    # end

def remove_wat(atomgroup):
    answer = brdige.AtomGroup()
    for grpkey, grp in atomgroup.groups():
        print(grp.name)
        if grpkey == 'HOH':
            print("remove HOH")
            continue
        elif grpkey == 'WAT':
            print("remove WAT")
            continue

        tmpgrp = remove_wat(grp)
        if tmpgrp.get_number_of_all_atoms() > 0:
            answer.set_group(grpkey, tmpgrp)
    for atmkey, atm in atomgroup.atoms():
        answer.set_atom(atmkey, atm)

    return answer



if __name__ == '__main__':
    main()
