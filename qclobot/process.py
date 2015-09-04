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

    def __init__(self, logfile_path = ''):
        self._procs = []
        self._logfile_path = str(logfile_path)

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
            
    def commit(self):
        stdout = self._procs[-1].stdout
        stderr = self._procs[-1].stderr
        stdout_fd = stdout.fileno()
        stderr_fd = stderr.fileno()
        self._make_non_blocking(stdout_fd)
        self._make_non_blocking(stderr_fd)

        #stdout_data = []
        #stderr_data = []
        stdout_eof = False
        stderr_eof = False

        logfile = None
        if len(self._logfile_path) > 0:
            logfile = open(self._logfile_path, mode='ab')
        
        while True:
            to_check = [stdout_fd]*(not stdout_eof) + [stderr_fd]*(not stderr_eof)
            ready = select.select(to_check, [], [])
            if stdout_fd in ready[0]:
                stdout_chunk = stdout.read()
                if stdout_chunk == '':
                    stdout_eof = True
                else:
                    sys.stdout.write(stdout_chunk)
                    sys.stdout.flush()
                    if logfile != None:
                        logfile.write(stdout_chunk)
                    #stdout_data.append(stdout_chunk)
            if stderr_fd in ready[0]:
                stderr_chunk = stderr.read()
                if stderr_chunk == '':
                    stderr_eof = True
                else:
                    sys.stderr.write(stderr_chunk)
                    sys.stderr.flush()
                    if logfile != None:
                        logfile.write(stdout_chunk)
                    #stderr_data.append(stderr_chunk)
            if stdout_eof and stderr_eof:
                break
            select.select([], [], [], .1)
            
        self._procs[-1].wait()
        status = self._procs[-1].returncode

        for p in self._procs:
            p.stdout.close()
        self._procs = []

        if logfile != None:
            logfile.close()
        
        return status


    def _make_non_blocking(self, fd):
        fl = fcntl.fcntl(fd, fcntl.F_GETFL)
        try:
            fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NDELAY)
            fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)
        except AttributeError:
            fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.FNDELAY)

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
