#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2002-2015 The ProteinDF project
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

import jinja2
import pprint

import proteindf_bridge as bridge

from .__version__ import __version__
from .qcframe import QcFrame
from .qcfragment import QcFragment
from .qccontrolobject import QcControlObject
from .qcerror import QcControlError


import logging

logger = logging.getLogger(__name__)


class QcControl(QcControlObject):
    _modeling = bridge.Modeling()

    def __init__(self):
        super(QcControl, self).__init__()
        self._cache = {}

        self._frames = {}
        self._frames["default"] = {}
        self._frames["default"]["guess"] = {"method": "harris", "force": False}
        self._frames["default"]["basis_set"] = "DZVP2"
        self._frames["default"]["brd_file"] = ""

        self._last_frame_name = ""

    # ------------------------------------------------------------------
    # property
    # ------------------------------------------------------------------
    @property
    def last_frame_name(self):
        """return the last frame name"""
        return self._last_frame_name

    def get_frame(self, name):
        return self._frames.get(name, None)

    # ------------------------------------------------------------------
    # task or frame
    # ------------------------------------------------------------------
    def _run_task_cmd(self, frame_data):
        assert isinstance(frame_data, dict)

        # condition ----------------------------------------------------
        if self._exec_condition_with_items(frame_data):
            return
        if self._exec_condition_when(frame_data):
            return

        # include ------------------------------------------------------
        retval_include_tasks = self._exec_include_tasks(frame_data)
        if retval_include_tasks:
            return retval_include_tasks

        # task ---------------------------------------------------------
        if self._exec_task_debug(frame_data):
            return

        if self._exec_task_mail(frame_data):
            return

        # frame --------------------------------------------------------
        if self._exec_frame_default(frame_data):
            return

        self._exec_frame_object(frame_data)

    # ------------------------------------------------------------------
    # condition
    # ------------------------------------------------------------------
    def _exec_condition_with_items(self, frame_data):
        is_break = False

        with_items = frame_data.get("with_items", None)
        if with_items:
            iter_items = list(self._vars.get(bridge.StrUtils.to_unicode(with_items), []))

            new_frame_data = dict(frame_data)
            new_frame_data.pop("with_items")
            yaml_str = bridge.get_yaml(new_frame_data)
            template = jinja2.Template(yaml_str)

            for item in iter_items:
                logger.info("template render: item={}".format(repr(item)))
                yaml_str = template.render(item=item)
                for new_frame_data in bridge.parse_yaml(yaml_str):
                    self._run_task(new_frame_data)
            is_break = True

        return is_break

    def _exec_condition_when(self, task_data):
        is_break = False

        when_phrase = task_data.get("when", None)
        if when_phrase:
            globals_data = self._vars
            judge = eval(when_phrase, globals_data)
            if judge:
                new_task_data = dict(task_data)
                new_task_data.pop("when")
                self._run_task(new_task_data)
            is_break = True

        return is_break

    # ------------------------------------------------------------------
    # task
    # ------------------------------------------------------------------
    def _exec_task_mail(self, in_task_data):
        is_break = False

        if "mail" in in_task_data:
            mail_data = in_task_data["mail"]
            mailer = bridge.Mail()
            mailer.smtp_server = mail_data.get("smtp_server")
            if "smtp_port" in mail_data:
                mailer.smtp_port = mail_data.get("smtp_port")
            if "use_SSL" in mail_data:
                mailer.use_SSL = mail_data.get("use_SSL")
            mailer.smtp_account = mail_data.get("smtp_account")
            mailer.smtp_password = mail_data.get("smtp_password")
            mailer.from_address = mail_data.get("from_address")
            mailer.to_address = mail_data.get("to_address")
            mailer.subject = mail_data.get("subject")
            mailer.text = mail_data.get("msg")
            mailer.send()

            is_break = True

        return is_break

    # ------------------------------------------------------------------
    # frame
    # ------------------------------------------------------------------
    def _exec_frame_default(self, in_frame_data):
        assert isinstance(in_frame_data, dict)
        is_break = False

        name = in_frame_data.get("name", None)
        if (name != None) and (name == "default"):
            default_data = dict(in_frame_data)
            default_data.pop("name")
            self._frames["default"].update(default_data)
            is_break = True

        return is_break

    def _exec_frame_object(self, frame_data):
        """
        execute something to frame object
        """
        assert isinstance(frame_data, dict)

        # --------------------------------------------------------------
        # frame name operation
        # --------------------------------------------------------------
        if "name" not in frame_data:
            logger.critical("name is not defined.")
        frame_name = frame_data.get("name")
        logger.info("#" * 72)
        logger.info("# FRAME: {} ".format(frame_name))
        logger.info("#" * 72)
        frame = QcFrame(frame_name)

        # --------------------------------------------------------------
        # setup default value
        # --------------------------------------------------------------
        frame_data = dict(self._frames["default"], **frame_data)
        # self._set_default(self._frames['default'], frame_data)

        # --------------------------------------------------------------
        # preparation of input
        # --------------------------------------------------------------
        # Lower compatibility: guess
        if isinstance(frame_data.get("guess", None), str):
            logger.debug("string-type is deprecated in 'guess' section. Use map(dict)-type.")
            guess_str = frame_data["guess"]
            frame_data["guess"] = {"method": guess_str, "force": False}

        # --------------------------------------------------------------
        # frame configuration
        # --------------------------------------------------------------
        if isinstance(frame_data["guess"], dict) != True:
            raise QcControlError(frame_data["guess"], "guess")
        assert isinstance(frame_data["guess"]["method"], str)
        frame.pdfparam.guess = frame_data["guess"]["method"]

        charge = frame_data.get("charge", None)
        if charge:
            frame.charge = charge

        frame.pdfparam.CDAM_tau = frame_data.get("CDAM_tau", None)
        frame.pdfparam.CD_epsilon = frame_data.get("CD_epsilon", None)

        frame.pdfparam.convergence_threshold_energy = frame_data.get("convergence/threshold_energy", None)
        frame.pdfparam.convergence_threshold = frame_data.get("convergence/threshold", None)
        frame.pdfparam.convergence_type = frame_data.get("convergence/type", None)

        frame.pdfparam.orbital_independence_threshold = frame_data.get("orbital_independence_threshold", None)
        frame.pdfparam.orbital_independence_threshold_canonical = frame_data.get(
            "orbital_independence_threshold/canonical", None
        )
        frame.pdfparam.orbital_independence_threshold_lowdin = frame_data.get(
            "orbital_independence_threshold/lowdin", None
        )

        frame.pdfparam.scf_acceleration = frame_data.get("scf_acceleration", None)
        frame.pdfparam.scf_acceleration_damping_start = frame_data.get(
            "scf_acceleration/damping/start", None
        )
        frame.pdfparam.scf_acceleration_damping_damping_factor = frame_data.get(
            "scf_acceleration/damping/damping_factor", None
        )
        frame.pdfparam.scf_acceleration_anderson_start = frame_data.get(
            "scf_acceleration/anderson/start", None
        )
        frame.pdfparam.scf_acceleration_anderson_damping_factor = frame_data.get(
            "scf_acceleration/anderson/damping_factor", None
        )

        frame.pdfparam.xc_functional = frame_data.get("XC_functional", None)
        frame.pdfparam.j_engine = frame_data.get("J_engine", None)
        frame.pdfparam.k_engine = frame_data.get("K_engine", None)
        frame.pdfparam.xc_engine = frame_data.get("XC_engine", None)

        frame.pdfparam.gridfree_dedicated_basis = frame_data.get("gridfree/dedicated_basis", None)
        frame.pdfparam.gridfree_orthogonalize_method = frame_data.get("gridfree/orthogonalize_method", None)
        frame.pdfparam.gridfree_CDAM_tau = frame_data.get("gridfree/CDAM_tau", None)
        frame.pdfparam.gridfree_CD_epsilon = frame_data.get("gridfree/CD_epsilon", None)

        frame.pdfparam.extra_keywords = frame_data.get("pdf_extra_keywords", {})
        # logger.info("pdf extra keywords:")
        # logger.info(repr(frame_data["pdf_extra_keywords"]))

        # properties -----------
        # cmd alias
        if "cmd_alias" in frame_data:
            frame.set_command_alias(frame_data.get("cmd_alias"))

        if "db_filename" in frame_data:
            frame.set_db_filename(frame_data.get("db_filename"))

        # fragments ------------
        logger.info("> make fragments for [{frame_name}]".format(frame_name=frame_name))
        fragments_list = self._get_fragments(frame_data.get("fragments", []), frame_data)
        for f in fragments_list:
            assert isinstance(f, QcFragment)
            frame[f.name] = f

            # fragment information
            # parent_frame_name = ""
            # if f.parent != None:
            #    parent_frame_name = f.parent.name
            # logger.info("> add fragment {ref_frame}/{ref_fragment} to {frame_name}".format(
            #    ref_frame=parent_frame_name,
            #    ref_fragment=f.name,
            #    frame_name=frame_name
            # ))
            # for subgrp_name, subgrp in frame[f.name].groups():
            #    logger.info('  > subfragment: {}, parent={}'.format(subgrp_name, subgrp.parent.name))
        logger.debug("> make fragments for [{frame_name}]: end".format(frame_name=frame_name))

        # --------------------------------------------------------------
        # action
        # --------------------------------------------------------------
        logger.info("::action")
        frame.save()
        self._frames[frame_name] = frame
        self._last_frame_name = frame_name

        # pre-SCF
        if frame_data.get("pre_scf", False):
            frame.calc_preSCF()

        # make guess by QCLObot
        logger.info("::GUESS")
        guess = frame_data.get("guess", {"method": "harris"})
        assert isinstance(guess, dict)

        guess_method = guess.get("method", "harris")
        guess_method = guess_method.lower()
        guess_force = guess.get("force", False)
        if isinstance(guess_force, str):
            guess_force = guess_force.lower()
            if (len(guess_force) > 0) and (guess_force in ("true", "yes", "1")):
                guess_force = True
            else:
                guess_force = False

        if guess_method == "density":
            frame.guess_density(force=guess_force)
        elif guess_method == "qclo":
            frame.guess_QCLO(force=guess_force)
        else:
            pass
            # frame.guess_Harris()

        # Single point calc.
        if frame_data.get("sp", False):
            frame.calc_sp()

        # force
        if frame_data.get("force", False):
            frame.calc_force()

        # summary
        summary_act = frame_data.get("summary", False)
        if isinstance(summary_act, dict):
            format_str = summary_act.get("format", None)
            filepath = summary_act.get("filepath", None)
            frame.summary(format_str=format_str, filepath=filepath)
        elif isinstance(summary_act, bridge.basestring):
            frame.summary(format_str=summary_act)
        elif isinstance(summary_act, bool):
            if summary_act:
                frame.summary()

        # finalize
        frame.save()

    # ------------------------------------------------------------------
    # utils
    # ------------------------------------------------------------------
    def _get_value(self, keyword, input_data):
        """
        return a value corresponding to the keyword from the input_data or defaults
        """
        keyword = bridge.StrUtils.to_unicode(keyword)

        answer = None
        if keyword in input_data:
            answer = input_data.get(keyword)
        elif keyword in self._frames["default"]:
            answer = self._frames["default"][keyword]
        return answer

    # ------------------------------------------------------------------
    # fragment
    # ------------------------------------------------------------------
    def _get_fragments(self, fragments_data, default, num_of_nest = 0):
        assert isinstance(fragments_data, list)
        indent = "  " * num_of_nest 
        logger.info("{}>>> make list of fragments: start".format(indent))
        logger.debug(str(default))

        # fragment = QcFragment()
        answer = []
        for frg_data in fragments_data:
            assert isinstance(frg_data, dict)

            self._set_default(default, frg_data)
            logger.debug("> get_fragments(): loop")
            logger.debug(frg_data)
            name = frg_data.get("name", "")
            logger.debug("> name: {}".format(name))

            subfrg = None
            if "fragments" in frg_data:
                logger.info("{indent}>>>> collect fragments for {name}: start".format(indent=indent, name=name))
                subfrg_list = self._get_fragments(frg_data.get("fragments"), default, num_of_nest +1)
                subfrg = QcFragment(name=name)
                for item in subfrg_list:
                    logger.debug("> list name: {}".format(item.name))
                    subfrg[item.name] = item
                logger.info("{indent}<<<< collect fragments for {name}: end".format(indent=indent, name=name))
            elif "add_H" in frg_data:
                subfrg = self._get_add_H(frg_data)
            elif "add_CH3" in frg_data:
                subfrg = self._get_add_CH3(frg_data)
            elif "add_ACE" in frg_data:
                subfrg = self._get_add_ACE(frg_data)
            elif "add_NME" in frg_data:
                subfrg = self._get_add_NME(frg_data)
            elif "brd_select" in frg_data:
                subfrg = self._get_default_fragment(frg_data)
            elif "atomlist" in frg_data:
                subfrg = self._getfrg_atomlist(frg_data)
            elif "reference" in frg_data:
                subfrg = self._get_reference_fragment(frg_data)
            else:
                logger.critical("unknown fragment command/object: {}".format(str(frg_data.keys())))
                raise QcControlError("fragment", "unknown fragment: {}".format(str(frg_data)))

            if subfrg == None:
                frg_data_str = pprint.pformat(frg_data)
                logger.critical("> unknown sub-fragment:")
                logger.critical(">  sub-fragment information:")
                logger.critical(frg_data_str)
                raise QcControlError("unknown subfragment:", frg_data_str)

            assert isinstance(subfrg, QcFragment)
            subfrg_atomgroup = subfrg.get_AtomGroup()
            logger.info(
                "{indent}subfrg append: {name} ({formula}: charge={charge})".format(
                    indent=indent,
                    name=name,
                    formula=subfrg_atomgroup.get_formula(),
                    charge=subfrg_atomgroup.charge,
                )
            )
            answer.append(subfrg)

        logger.info("{indent}<<<< make list of fragments: end".format(indent=indent))
        return answer

    def _set_default(self, default_values, update_values):
        """
        copy default values to variable
        """
        keywords = [
            "brd_file",
            "basis_set",
            "basis_set_aux",
            "basis_set_gridfree",
            "guess",
            "cut_value",
            "CDAM_tau",
            "CD_epsilon",
            "orbital_independence_threshold",
            "orbital_independence_threshold/canonical",
            "orbital_independence_threshold/lowdin",
            "gridfree/dedicated_basis",
            "gridfree/orthogonalize_method",
            "gridfree/CDAM_tau",
            "gridfree/CD_epsilon",
            "scf_acceleration",
            "scf_acceleration/damping/damping_factor",
            "scf_acceleration/anderson/start",
            "scf_acceleration/anderson/damping_factor",
            "convergence/threshold_energy",
            "convergence/threshold",
            "XC_functional",
            "J_engine",
            "K_engine",
            "XC_engine",
            "pdf_extra_keywords",
        ]
        for keyword in keywords:
            update_values.setdefault(keyword, default_values.get(keyword, None))

    def _get_default_fragment(self, frg_data):
        assert isinstance(frg_data, dict)

        brd_file_path = frg_data.get("brd_file")

        brd_select = frg_data.get("brd_select")

        select_path = ""
        except_path = ""
        if isinstance(brd_select, str):
            select_path = brd_select
        elif isinstance(brd_select, dict):
            # key: path
            if "path" in brd_select:
                select_path = brd_select.get("path", "")
            else:
                raise QcControlError(frg_data, "not found path key in brd_select.")

            # key: except
            except_path = brd_select.get("except", "")

        else:
            raise QcControlError(frg_data, "unknown brd_select type.")

        # select atomgroup
        atomgroup = self._select_atomgroup(brd_file_path, select_path)
        if len(except_path) > 0:
            except_atomgroup = self._select_atomgroup(brd_file_path, except_path)

            except_atomgroup = atomgroup & except_atomgroup
            atomgroup = atomgroup ^ except_atomgroup

        frg = QcFragment(atomgroup)
        frg.margin = False
        if "name" not in frg_data:
            raise QcControlError("name keyword NOT FOUND", str(frg_data))
        frg.name = frg_data.get("name")

        self._set_basis_set(
            frg,
            frg_data.get("basis_set"),
            frg_data.get("basis_set_aux", None),
            frg_data.get("basis_set_gridfree", None),
        )

        return frg

    def _get_reference_fragment(self, frg_data):
        """return new QcFragment object referenced from input(frg_data)"""
        assert isinstance(frg_data, dict)

        if not "frame" in frg_data["reference"]:
            raise QcControlError("NOT FOUND frame key in reference fragment", pprint.pformat(frg_data))
        if not "fragment" in frg_data["reference"]:
            raise QcControlError("NOT FOUND fragment key in reference fragment", pprint.pformat(frg_data))

        ref_frame_name = bridge.StrUtils.to_unicode(frg_data["reference"]["frame"])
        ref_fragment_name = bridge.StrUtils.to_unicode(frg_data["reference"]["fragment"])

        if ref_frame_name not in self._frames:
            raise QcControlError("UNKNOWN FRAME", ref_frame_name)
        ref_frame = self._frames.get(ref_frame_name)

        if not ref_frame.has_fragment(ref_fragment_name):
            logger.critical(str(self._frames[ref_frame_name]))
            raise QcControlError("UNKNOWN FRAGMENT", "{} in {}".format(ref_fragment_name, ref_frame_name))

        logger.info("> reference fragment: [{}/{}]".format(ref_frame_name, ref_fragment_name))
        ref_fragment = ref_frame[ref_fragment_name]

        fragment = QcFragment(ref_fragment)
        fragment.name = frg_data.get("name")
        fragment.ref_fragment = ref_fragment
        return fragment

    def _get_add_H(self, frg_data):
        assert isinstance(frg_data, dict)
        brd_file_path = frg_data.get("brd_file")
        brd_select_H = frg_data.get("displacement")
        atomgroup_H = self._select_atomgroup(brd_file_path, brd_select_H)
        atomgroup_H = atomgroup_H.get_atom_list()
        assert len(atomgroup_H) > 0
        atom_H = atomgroup_H[0]
        atom_H.symbol = "H"

        ag_H = bridge.AtomGroup()
        ag_H.set_atom("H*", atom_H)

        H = QcFragment(ag_H)
        H.margin = True
        if "name" not in frg_data:
            raise QcControlError("NOT FOUND name key in add_H", str(frg_data))
        H.name = frg_data.get("name")

        self._set_basis_set(
            H,
            frg_data.get("basis_set"),
            frg_data.get("basis_set_aux", None),
            frg_data.get("basis_set_gridfree", None),
        )

        return H

    def _get_add_CH3(self, frg_data):
        assert isinstance(frg_data, dict)
        brd_file_path = frg_data.get("brd_file")
        brd_select_C1 = frg_data.get("displacement")
        brd_select_C2 = frg_data.get("root")
        logger.info("add CH3: displacement: {}, root: {}".format(brd_select_C1, brd_select_C2))

        atomgroup_C1 = self._select_atomgroup(brd_file_path, brd_select_C1)
        atomgroup_C1 = atomgroup_C1.get_atom_list()
        if len(atomgroup_C1) == 0:
            logger.critical("cannnot find displacement atom: key={}".format(brd_select_C1))
            raise QcControlError('No atoms found for "displacement" in add_CH3', brd_select_C1)
        atom_C1 = atomgroup_C1[0]

        atomgroup_C2 = self._select_atomgroup(brd_file_path, brd_select_C2)
        atomgroup_C2 = atomgroup_C2.get_atom_list()
        if len(atomgroup_C2) == 0:
            logger.critical("cannnot find displacement atom: key={}".format(brd_select_C2))
            raise QcControlError('No atoms found for "root" in add_CH3')
        atom_C2 = atomgroup_C2[0]
        ag_CH3 = self._modeling.add_methyl(atom_C1, atom_C2)
        ag_CH3.set_atom("C", atom_C1)

        CH3 = QcFragment(ag_CH3)
        CH3.margin = True
        if "name" not in frg_data:
            raise QcControlError("NOT FOUND name key in add_CH3", str(frg_data))
        CH3.name = frg_data.get("name")

        self._set_basis_set(CH3, frg_data.get("basis_set"), frg_data.get("basis_set_gridfree", None))

        return CH3

    def _get_add_ACE(self, frg_data):
        assert isinstance(frg_data, dict)
        brd_file_path = frg_data.get("brd_file")
        brd_select = frg_data.get("brd_select")
        atomgroup = self._select_atomgroup(brd_file_path, brd_select)
        ag_ACE = self._modeling.get_ACE_simple(atomgroup)

        ACE = QcFragment(ag_ACE)
        ACE.margin = True
        if "name" not in frg_data:
            raise QcControlError("NOT FOUND name key in add_ACE", str(frg_data))
        ACE.name = frg_data.get("name")

        self._set_basis_set(
            ACE,
            frg_data.get("basis_set"),
            frg_data.get("basis_set_aux", None),
            frg_data.get("basis_set_gridfree", None),
        )

        return ACE

    def _get_add_NME(self, frg_data):
        assert isinstance(frg_data, dict)
        brd_file_path = frg_data.get("brd_file")
        brd_select = frg_data.get("brd_select")
        atomgroup = self._select_atomgroup(brd_file_path, brd_select)
        ag_NME = self._modeling.get_NME_simple(atomgroup)

        NME = QcFragment(ag_NME)
        NME.margin = True
        if "name" not in frg_data:
            raise QcControlError("NOT FOUND name key in add_NME", str(frg_data))
        NME.name = frg_data.get("name")

        self._set_basis_set(
            NME,
            frg_data.get("basis_set"),
            frg_data.get("basis_set_aux", None),
            frg_data.get("basis_set_gridfree", None),
        )

        return NME

    def _getfrg_atomlist(self, frg_data):
        assert isinstance(frg_data, dict)
        atomlist = frg_data.get("atomlist")

        atomgroup = bridge.AtomGroup()
        if "name" not in frg_data:
            raise QcControlError("NOT FOUND name key in atomlist", str(frg_data))
        atomgroup.name = frg_data.get("name")

        index = 1
        if isinstance(atomlist, list):
            for line in atomlist:
                if isinstance(line, bridge.basestring):
                    if line[0] == "/":
                        # case '/model/...'
                        brd_select = line
                        brd_file_path = frg_data.get("brd_file")
                        selected_atomgroup = self._select_atomgroup(brd_file_path, brd_select)
                        # logger.debug('select path: {}'.format(brd_select))
                        # logger.debug('selected: {}'.format(str(selected_atomgroup)))
                        if selected_atomgroup.get_number_of_all_atoms() > 0:
                            for selected_atom in selected_atomgroup.get_atom_list():
                                atomgroup.set_atom(index, selected_atom)
                                index += 1
                        else:
                            logger.warn("not math: {}".format(brd_select))
                        continue
                    else:
                        # case: "H 1.00 2.00 3.00"
                        line = line.split()
                if isinstance(line, list):
                    if len(line) == 4:
                        # case [H, 1.00, 2.00, 3.00]
                        atom = bridge.Atom()
                        atom.symbol = line[0]
                        atom.xyz = bridge.Position(line[1], line[2], line[3])
                        atomgroup.set_atom(index, atom)
                        index += 1
                    elif len(line) == 5:
                        # case [X, 1,00, 2.00, 3.00, -1.0] for dummy charge
                        atom = bridge.Atom()
                        atom.symbol = line[0]
                        atom.xyz = bridge.Position(line[1], line[2], line[3])
                        atom.charge = line[4]
                        atomgroup.set_atom(index, atom)
                        index += 1
                    else:
                        logger.error("mismatch atom data size: {}".format(str(line)))
                else:
                    logger.error("atomlist item shuld be array or string: {}".format(str(line)))
        else:
            logger.error("atomlist shuld be list: {}".format(str(atomlist)))

        frg = QcFragment(atomgroup)
        self._set_basis_set(
            frg,
            frg_data.get("basis_set"),
            frg_data.get("basis_set_aux", None),
            frg_data.get("basis_set_gridfree", None),
        )
        return frg

    def _select_atomgroup(self, brd_file_path, brd_select):
        brd_file_path = bridge.StrUtils.to_unicode(brd_file_path)
        brd_select = bridge.StrUtils.to_unicode(brd_select)
        # logger.debug('_get_atomgroup(): brd_path={}, select={}'.format(brd_file_path, brd_select))

        self._cache.setdefault("brdfile", {})
        if brd_file_path not in self._cache["brdfile"]:
            atomgroup = bridge.load_atomgroup(brd_file_path)

            self._cache["brdfile"][brd_file_path] = {}
            self._cache["brdfile"][brd_file_path]["atomgroup"] = atomgroup

        atomgroup = self._cache["brdfile"][brd_file_path]["atomgroup"]
        # selecter = bridge.Select_PathRegex(brd_select)
        selecter = bridge.Select_Path_wildcard(brd_select)
        answer = atomgroup.select(selecter)

        return answer

    def _set_basis_set(self, fragment, basis_set, basis_set_aux=None, basis_set_gridfree=None):
        """
        set basis set to the fragment object.
        """
        if isinstance(fragment, QcFragment):
            for subfrg_name, subfrg in fragment.groups():
                self._set_basis_set(subfrg, basis_set, basis_set_aux, basis_set_gridfree)
            for atm_name, atm in fragment.atoms():
                atm.basisset = "O-{}.{}".format(self._find_basis_set_name(atm, basis_set), atm.symbol)
                if basis_set_aux:
                    atm.basisset_j = "A-{}.{}".format(self._find_basis_set_name(atm, basis_set_aux), atm.symbol)
                    atm.basisset_xc = "A-{}.{}".format(self._find_basis_set_name(atm, basis_set_aux), atm.symbol)
                else:
                    atm.basisset_j = "A-{}.{}".format(self._find_basis_set_name(atm, basis_set), atm.symbol)
                    atm.basisset_xc = "A-{}.{}".format(self._find_basis_set_name(atm, basis_set), atm.symbol)
                if basis_set_gridfree:
                    atm.basisset_gridfree = "O-{}.{}".format(
                        self._find_basis_set_name(atm, basis_set_gridfree), atm.symbol
                    )
                else:
                    atm.basisset_gridfree = "O-{}.{}".format(self._find_basis_set_name(atm, basis_set), atm.symbol)
        else:
            raise QcControlError("Program Error: ", str(fragment))

    def _find_basis_set_name(self, atom, input_obj):
        """
        find suitable basis set name for the atom.

        """
        assert isinstance(atom, bridge.Atom)
        ans = ""
        if isinstance(input_obj, bridge.basestring):
            ans = bridge.StrUtils.to_unicode(input_obj)
        elif isinstance(input_obj, dict):
            if atom.name in input_obj:
                ans = input_obj[atom.name]
            elif atom.symbol in input_obj:
                ans = input_obj[atom.symbol]
            else:
                raise QcControlError("type mismatch.", str(input_obj))
        return ans

    # ------------------------------------------------------------------
    # others
    # ------------------------------------------------------------------

    def show_version(self):
        logger.info("=" * 80)
        logger.info("QCLObot version: {version}".format(version=str(__version__)))
        logger.info("=" * 80)


if __name__ == "__main__":
    pass
