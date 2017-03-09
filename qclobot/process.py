#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014-2015 The ProteinDF development team.
# see also AUTHORS and README if provided.
# 
# This file is a part of the ProteinDF software package.
# 
# The ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# The ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import io
import shlex
import subprocess
import fcntl
import select

import logging
logger = logging.getLogger(__name__)

import pdfbridge

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

    def cmd(self, cmd):
        p_args = {
            'stdin': None,
            'stdout': subprocess.PIPE,
            'stderr': subprocess.PIPE,
        }
        p_args['shell'] = self._is_shell_script(cmd)
        
        return self._exec(cmd, p_args)

    def pipe(self, cmd):
        p_args = {
            'stdout': subprocess.PIPE,
            'stderr': subprocess.PIPE,
        }
        p_args['shell'] = self._is_shell_script(cmd)

        if len(self._procs) == 0:
            p_args['stdin'] = None
        else:
            p_args['stdin'] = self._procs[-1].stdout

        return self._exec(cmd, p_args)

    def _exec(self, cmd, p_args):
        p_args['bufsize'] = -1

        if p_args.get('shell', False) == False:
            cmd = shlex.split(cmd)
            
        new_proc = None
        try:
            new_proc = subprocess.Popen(cmd, **p_args)
        except OSError as e:
            sys.stderr.write('Failed to execute command: {}\n'.format(cmd))
            raise
            # raise e
        except:
            raise

        self._procs.append(new_proc)
        return self
            
    def commit(self,
               stdout_filepath=None, stderr_filepath=None,
               stdout_through=True, stderr_through=True):
        stdout = self._procs[-1].stdout
        stderr = self._procs[-1].stderr

        stdout_file = None
        if stdout_filepath != None:
            stdout_file = open(stdout_filepath, mode='a')
        stderr_file = None
        if stderr_filepath != None:
            stderr_file = open(stderr_filepath, mode='a')

        # see also.
        # http://qiita.com/FGtatsuro/items/0f68ab9c1bcad9c4b320
        # https://gist.github.com/mattbornski/3299031
        with io.open(stdout.fileno(), closefd=False) as out_stream, io.open(stderr.fileno(), closefd=False) as err_stream:
            for line in out_stream:
                #line = line.rstrip('\n')
                #line = pdfbridge.Utils.to_unicode(line)
                if stdout_through:
                    sys.stdout.write(line)
                if stdout_file != None:
                    stdout_file.write(line)
            for line in err_stream:
                #line = line.rstrip('\n')
                #line = pdfbridge.Utils.to_unicode(line)
                if stderr_through:
                    sys.stderr.write(line)
                if stderr_file != None:
                    stderr_file.write(line)

            
        self._procs[-1].wait()
        status = self._procs[-1].returncode

        for p in self._procs:
            p.stdout.close()
            p.stderr.close()
        self._procs = []

        if stdout_file != None:
            stdout_file.close()
        if stderr_file != None:
            stderr_file.close()
            
        
        return status

    def _is_shell_script(self, cmd):
        answer = False

        cmd_list = shlex.split(cmd)
        if len(cmd_list) > 0:
            name, ext = os.path.splitext(cmd_list[0])
            if ext == '.sh':
                answer = True

        return answer
            
            
if __name__ == '__main__':
    import doctest
    doctest.testmod()
