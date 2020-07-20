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


class QcError(Exception):
    """
    Base class for QCLO module
    """

    def __init__(self, errmsg=""):
        self.errmsg = errmsg

    def __str__(self):
        return "QCLObot error: {}".format(self.errmsg)


class QcControlError(QcError):
    """
    Exception raised for errors in the control input.

    Attributes:
        expr (str): input expression in which the error occurred
        msg (str): explanation of the error
    """

    def __init__(self, expr, msg):
        self.errmsg = "Input Error: {} ({})".format(msg, str(expr))


class QcScriptRunningError(QcError):
    """
    Exception raised for errors to run script.

    Attributes:
        cmd_name (str, list): command name
        err_msg (str): error message on running
    """

    def __init__(self, cmd_name, err_msg):
        if isinstance(cmd_name, list):
            cmd_name = ' '.join(cmd_name)
        self.errmsg = 'The following script could not run: "{}"\nerr_msg:{}'.format(
            cmd_name, err_msg)


class QcTaskError(QcError):
    def __inti__(self, expr, msg):
        self.errmsg = "Input Error: {} ({})".format(msg, str(expr))
