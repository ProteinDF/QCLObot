#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import io
import shlex
import subprocess
import fcntl
import select

import proteindf_bridge as bridge

import logging
logger = logging.getLogger(__name__)


class Process(object):
    '''
    create Process object.
    >>> p = Process()

    register command.
    >>> p.cmd('echo "Hello world!"')
    ... #doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    <__main__.Process object at 0x...>

    use pipe.
    >>> p.pipe('sed -e "s/world/python/"')
    ... #doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    <__main__.Process object at 0x...>

    execute command.
    >>> return_code = p.commit()
    ... #doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    Hello python!

    ...
    '''

    def __init__(self):
        self._procs = []
        # logger.info("DEFAULT_BUFFER_SIZE={}".format(io.DEFAULT_BUFFER_SIZE))

    def _is_shell_script(self, cmd):
        assert(isinstance(cmd, (str, list)))
        return isinstance(cmd, str)

    def run(self, cmd, stdout_filepath, stderr_filepath):
        p_args = {
            'stdin': None,
            'stdout': subprocess.PIPE,
            'stderr': subprocess.PIPE,
            'shell': self._is_shell_script(cmd)
        }
        p = subprocess.Popen(cmd, **p_args)
        (stdout, stderr) = p.communicate()

    def cmd(self, cmd):
        p_args = {
            'stdin': None,
            'stdout': subprocess.PIPE,
            'stderr': subprocess.PIPE,
            'shell': self._is_shell_script(cmd)
        }

        return self._exec(cmd, p_args)

    def pipe(self, cmd):
        p_args = {
            'stdout': subprocess.PIPE,
            'stderr': subprocess.PIPE,
            'shell': self._is_shell_script(cmd)
        }

        if len(self._procs) == 0:
            p_args['stdin'] = None
        else:
            p_args['stdin'] = self._procs[-1].stdout

        return self._exec(cmd, p_args)

    def _exec(self, cmd, p_args):
        # p_args['bufsize'] = -1  # use io.DEFAULT_BUFFER_SIZE (default)

        # if p_args['shell'] is False:
        # cmd = shlex.split(cmd)

        new_proc = None
        try:
            new_proc = subprocess.Popen(cmd, **p_args)
        except OSError as e:
            sys.stderr.write('Failed to execute command: {}\n'.format(str(cmd)))
            raise e

        self._procs.append(new_proc)
        return self

    def commit_blocked(self,
                       stdout_filepath=None, stderr_filepath=None,
                       stdout_through=True, stderr_through=True):
        (stdout_data, stderr_data) = self._procs[-1].communicate()
        status = self._procs[-1].returncode

        for p in self._procs:
            p.stdout.close()
            p.stderr.close()
        self._procs = []

        if stdout_filepath is not None:
            with open(stdout_filepath, mode='ab') as stdout_file:
                stdout_file.write(stdout_data)
        if stderr_filepath is not None:
            with open(stderr_filepath, mode='ab') as stderr_file:
                stderr_file.write(stderr_data)

        return status


if __name__ == '__main__':
    import doctest
    doctest.testmod()
