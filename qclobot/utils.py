#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)

import pdfbridge

def get_model(models):
    assert(check_format_models(models))

    model = None
    if models.get_number_of_groups() > 0:
        model_ids = models.get_group_list()
        model = models.get_group(model_ids.pop())
        
    else:
        logger.critical("not found any model.")
        raise
    return model


def check_format_models(models):
    assert(isinstance(models, pdfbridge.AtomGroup))
    answer = True

    if models.get_number_of_groups() > 0:
        for model_id, model in models.groups():
            answer &= check_format_model(model)
    else:
        logger.warning("no groups found in models: name={}".format(models.name))
    
    if models.get_number_of_atoms() != 0:
        logger.error("not allowed any atoms in models: name={}".format(models.name))
        answer = False
    return answer


def check_format_model(model):
    assert(isinstance(model, pdfbridge.AtomGroup))
    answer = True

    if model.get_number_of_groups() > 0:
        for chain_key, chain in model.groups():
            answer &= check_format_chain(chain)
    else:
        logger.warning("no groups found in model: name={}".format(model.name))

    if model.get_number_of_atoms() != 0:
        logger.error("not allowed any atoms in model: name={}".format(model.name))
        answer = False
    
    return answer


def check_format_chain(chain):
    assert(isinstance(chain, pdfbridge.AtomGroup))
    answer = True

    if chain.get_number_of_groups() > 0:
        for res_key, res in chain.groups():
            answer &= check_format_residue(res)
    else:
        logger.warning("no groups found in chain: name={}".format(chain.name))
    
    if chain.get_number_of_atoms() != 0:
        logger.error("not allowed any atoms in chain: name={}".format(chain.name))
        answer = False
        
    return answer
    

def check_format_residue(res):
    assert(isinstance(res, pdfbridge.AtomGroup))
    answer = True

    if res.get_number_of_groups() > 0:
        logger.error("not allowed groups in residue: name={}".format(res.name))
        answer = False

    if res.get_number_of_atoms() == 0:
        logger.warning("no atoms found in residue: name={}".format(res.name))
    
    return answer


def remove_WAT(atomgroup):
    """ remove water(WAT or HOH) residues
    """
    assert(isinstance(atomgroup, pdfbridge.AtomGroup))

    answer = pdfbridge.AtomGroup(atomgroup)
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
        

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    
