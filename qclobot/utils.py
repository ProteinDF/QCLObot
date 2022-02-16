#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import inspect

import proteindf_bridge as bridge

import logging

logger = logging.getLogger(__name__)


def locate(depth=0):
    """
    return tuple of (file, function, line number)
    cf.) https://qiita.com/ymko/items/b46d32b98f013f06d805
    """
    frame = inspect.currentframe().f_back
    return (
        os.path.basename(frame.f_code.co_filename),
        frame.f_code.co_name,
        frame.f_lineno,
    )


def file2atomgroup(input_path):
    assert isinstance(input_path, str)
    atomgroup = None

    abspath = os.path.abspath(input_path)
    (basename, ext) = os.path.splitext(abspath)
    ext = ext.lower()
    if ext in (".pdb", ".ent"):
        logger.info("load {path} as PDB text file.".format(path=abspath))
        pdb_obj = bridge.Pdb()
        pdb_obj.load(abspath)
        atomgroup = pdb_obj.get_atomgroup()
    else:
        logger.info("load {path} as bridge file.".format(path=abspath))
        atomgroup = bridge.load_atomgroup(abspath)

    return atomgroup


def get_model(models):
    assert bridge.Format.is_models(models)

    model = None
    if models.get_number_of_groups() > 0:
        model_ids = models.get_group_list()
        model = models.get_group(model_ids.pop())

    else:
        logger.critical("not found any model.")
        raise
    return model


def find_max_chain_id(model):
    """return max chain ID"""
    answer = ord("A")
    if bridge.Format.is_protein(model):
        for chain_id, chain in model.groups():
            answer = max(answer, ord(chain_id))
    return chr(answer)


# TODO: use bridge.Utils.remove_WAT()
def remove_WAT(atomgroup):
    """remove water(WAT or HOH) residues"""
    assert isinstance(atomgroup, bridge.AtomGroup)

    answer = bridge.AtomGroup(atomgroup)
    wat_keys = ["HOH", "WAT"]
    remove_groups = []
    for key, grp in answer.groups():
        grp_name = grp.name
        logger.debug("check group name: {}".format(grp_name))
        if grp_name in wat_keys:
            logger.debug("remove name: {}".format(grp_name))
            remove_groups.append(key)
            continue

        answer.set_group(key, remove_WAT(grp))

    for key in remove_groups:
        answer.remove_group(key)

    return answer


if __name__ == "__main__":
    import doctest

    doctest.testmod()
