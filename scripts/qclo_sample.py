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

import argparse
import pprint

import proteindf_bridge as bridge
import proteindf_tools as pdf
import qclobot as qclo

import logging
import logging.config


def main():
    parser = argparse.ArgumentParser(description='QCLO program')
    parser.add_argument('brdfile',
                        nargs=1,
                        help='Bridge file')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    brdfile_path = args.brdfile[0]

    # load
    if verbose:
        print("reading: {}".format(brdfile_path))

    # manager = qclo.QcFrameManager()

    # all models
    models = bridge.load_atomgroup(brdfile_path)
    #print('>>>> models')
    # print(models)
    # print('----')

    # model
    model = models["model_1"]
    # print(model)

    frames = {}
    #
    N_TERM = 0  # N-termがNMEなら1
    C_TERM = 0  # C-termがACEなら1
    modeling = bridge.Modeling()
    for chain_name, chain in model.groups():
        max_resid = chain.get_number_of_groups()
        logging.info("chain {}: max_resid={}".format(chain_name, max_resid))
        # logging.info(str(chain))

        # for resid, res in chain.groups():
        for resid in range(1, 6):
            #resid = int(resid)
            res = chain[resid]

            logging.info(">>>> resid: {}".format(resid))
            res_model = bridge.AtomGroup(res)

            logging.debug("resid = {}".format(resid))

            # create base fragment
            fragment_res = qclo.QcFragment(res_model)

            ACE = qclo.QcFragment()
            if resid > 1 + N_TERM:
                #logging.debug("add ACE")

                prev = chain[resid - 1]
                #logging.debug(">>> prev")
                #logging.debug("\n" + str(prev))
                # ag_ACE = modeling.get_ACE(res_model, prev)
                ag_ACE = modeling.get_ACE_simple(prev)
                ACE = qclo.QcFragment(ag_ACE)

            NME = qclo.QcFragment()
            if resid + C_TERM < max_resid:
                #logging.debug("add NME")

                ahead = chain[resid + 1]
                #logging.debug(">>> ahead")
                #logging.debug("\n" + str(ahead))
                # ag_NME = modeling.get_NME(res_model, ahead)
                ag_NME = modeling.get_NME_simple(ahead)
                NME = qclo.QcFragment(ag_NME)

            name = 'res_{}-{}'.format(resid, resid)
            frame = qclo.QcFrame(name=name)
            frame[name] = fragment_res
            frame['ACE'] = ACE
            frame['NME'] = NME

            frames[frame.name] = frame
            setup_calc_conf(frame.pdfparam)
            frame.calc_sp(dry_run=False)
            frame.pickup_density_matrix()

        # step2
        # for resid, res in chain.groups():
        for resid in range(1, 6):
            resid = int(resid)
            if (resid + 2) > chain.get_number_of_groups():
                break
            calc3res(frames, chain, resid)

        # step3
        calc5res(frames, chain, 1)

        calc_test(frames)

    exit()


def calc3res(frames, chain, start):
    res1 = start
    res2 = start + 1
    res3 = start + 2

    if res1 >= 4:
        return

    frame = qclo.QcFrame(name='res_{}-{}'.format(res1, res3))

    res1str = 'res_{res1}-{res1}'.format(res1=res1)
    res2str = 'res_{res2}-{res2}'.format(res2=res2)
    res3str = 'res_{res3}-{res3}'.format(res3=res3)
    frame[res1str] = frames[res1str][res1str]
    frame[res2str] = frames[res2str][res2str]
    frame[res3str] = frames[res3str][res3str]

    frame['ACE'] = frames[res1str]['ACE']
    frame['NME'] = frames[res3str]['NME']

    frames[frame.name] = frame

    try:
        setup_calc_conf(frame.pdfparam)
        frame.guess_density()
        frame.calc_sp(dry_run=False)
        frame.calc_lo()
        frame.pickup_lo()
    except:
        raise


def calc5res(frames, chain, start):
    res1 = start
    res2 = start + 1
    res3 = start + 2
    res4 = start + 3
    res5 = start + 4

    frame = qclo.QcFrame(name='res_{}-{}'.format(1, 5))

    res1str = 'res_{res1}-{res1}'.format(res1=res1)
    res2str = 'res_{res2}-{res2}'.format(res2=res2)
    res3str = 'res_{res3}-{res3}'.format(res3=res3)
    res4str = 'res_{res4}-{res4}'.format(res4=res4)
    res5str = 'res_{res5}-{res5}'.format(res5=res5)

    res1_3str = 'res_{res1}-{res3}'.format(res1=res1, res3=res3)
    res2_4str = 'res_{res2}-{res4}'.format(res2=res2, res4=res4)
    res3_5str = 'res_{res3}-{res5}'.format(res3=res3, res5=res5)

    res1_5str = 'res_{res1}-{res5}'.format(res1=res1, res5=res5)

    frame[res1_5str] = qclo.QcFragment()
    frame[res1_5str].set_group('res1', frames[res1_3str][res1str])
    frame[res1_5str].set_group('res2', frames[res1_3str][res2str])
    frame[res1_5str].set_group('res3', frames[res2_4str][res3str])
    frame[res1_5str].set_group('res4', frames[res3_5str][res4str])
    frame[res1_5str].set_group('res5', frames[res3_5str][res5str])

    frame['ACE'] = frames['res_1-3']['ACE']
    frame['NME'] = frames['res_3-5']['NME']

    frames[frame.name] = frame

    try:
        setup_calc_conf(frame.pdfparam)
        frame.pdfparam.step_control = 'integral'
        setup_calc_conf(frame.pdfparam)
        frame.calc_preSCF(dry_run=False)
        frame.guess_QCLO()
        # frame.calc_sp()
    except:
        raise


def calc_test(frames):
    frame = qclo.QcFrame(name='res_1-3.QCLO')

    frame['res_1-3'] = qclo.QcFragment()
    frame['res_1-3'].set_group('res1', frames['res_1-3']['res_1-1'])
    frame['res_1-3'].set_group('res2', frames['res_1-3']['res_2-2'])
    frame['res_1-3'].set_group('res3', frames['res_1-3']['res_3-3'])

    frame['NME'] = frames['res_1-3']['NME']

    frames[frame.name] = frame

    try:
        setup_calc_conf(frame.pdfparam)
        frame.calc_preSCF()
        frame.guess_QCLO()
        # frame.calc_sp()
    except:
        raise


def setup_calc_conf(pdfparam):
    pdfparam.j_engine = 'CD'
    pdfparam.k_engine = 'CD'
    pdfparam.xc_engine = 'gridfree_CD'
    pdfparam.xc_functional = 'b3lyp'


if __name__ == '__main__':
    logging.config.fileConfig('logconfig.ini')
    main()
