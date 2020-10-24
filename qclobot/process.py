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

    def _set_nonblocking(self, fh):
        """ ファイルハンドルをnonblockingにする """
        fd = fh.fileno()
        fl = fcntl.fcntl(fd, fcntl.F_GETFL)
        fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)

    def _nonblocking_pass(self, infh, out_filehandles, blocking=False):
        """ nonblockingな入力ファイルハンドルinfhからデータを読み込んでそれをそのまま出力ファイルハンドルに渡す。
            catのようなパイプ。blocking==Trueだと普通のblockingのような振舞に。
        """
        assert(isinstance(out_filehandles, list))

        while True:
            try:
                buf = infh.read(4096)
                print(buf)
            except OSError as e:
                if e.args[0] == 11:   # Resource temporarily unavailable
                    if blocking:
                        # データが無いと再びトライ。selectでデータが来るのを待った方が無難かも。
                        continue
                else:
                    # nonblockingモード。
                    # 今はデータが無いので、他のIOができるようにループから出るようにする。
                    buf = ''

            if buf == '':
                # データが無いので、コントロールをcallerに返す。
                break

            for o_fh in out_filehandles:
                o_fh.write(buf)

    def cmd(self, cmd, stderr=None):
        """[summary]

        Args:
            cmd ([type]): [description]
            stderr (str): [description]. Defaults to None.

        Returns:
            [type]: [description]
        """
        # 'stderr': subprocess.PIPE, # stderrがPIPEに流れると刺さる。コンストラクタで作成されるとあとで変更できない
        if isinstance(stderr, str):
            stderr = stderr.upper()
            if stderr == "NULL":
                stderr = subprocess.DEVNULL
            elif stderr == "PIPE":
                stderr = subprocess.PIPE

        p_args = {
            'stdin': None,
            'stdout': subprocess.PIPE,
            'stderr': stderr,
            'shell': self._is_shell_script(cmd),
        }

        return self._exec(cmd, p_args)

    def pipe(self, cmd, stderr=None):
        """[summary]

        Args:
            cmd ([type]): [description]
            stderr ([type], optional): [description]. Defaults to None.

        Returns:
            [type]: [description]
        """
        if isinstance(stderr, str):
            stderr = stderr.upper()
            if stderr == "NULL":
                stderr = subprocess.DEVNULL
            elif stderr == "PIPE":
                stderr = subprocess.PIPE

        p_args = {
            'stdout': subprocess.PIPE,
            'stderr': stderr,
            'shell': self._is_shell_script(cmd),
        }

        if len(self._procs) == 0:
            p_args['stdin'] = None
        else:
            p_args['stdin'] = self._procs[-1].stdout

        return self._exec(cmd, p_args)

    def _exec(self, cmd, p_args):
        # p_args['bufsize'] = -1  # use io.DEFAULT_BUFFER_SIZE (default)
        # p_args['bufsize'] = 0  # no buffer
        # p_args['bufsize'] = 1  # line buffer
        # p_args['bufsize'] = 8192000

        # if p_args.get('shell', False) is False:
        #     cmd = shlex.split(cmd)

        new_proc = None
        try:
            new_proc = subprocess.Popen(cmd, **p_args)
        except OSError as e:
            cmd_str = cmd
            if isinstance(cmd, list):
                cmd_str = ' '.join(cmd)
            sys.stderr.write('Failed to execute command: {}\n'.format(cmd_str))
            raise e

        self._procs.append(new_proc)
        return self

    def commit_blocked(self,
                       stdout_filepath=None, stderr_filepath=None,
                       stdout_through=True, stderr_through=True):
        (stdout_data, stderr_data) = self._procs[-1].communicate()
        status = self._procs[-1].returncode

        # for p in self._procs:
        #     if p.stdout is not None:
        #         p.stdout.close()
        #     if p.stderr is not None:
        #         p.stderr.close()
        self._procs = []

        if stdout_filepath is not None:
            with open(stdout_filepath, mode='ab') as stdout_file:
                if stdout_data is not None:
                    stdout_file.write(stdout_data)
        if stderr_filepath is not None:
            with open(stderr_filepath, mode='ab') as stderr_file:
                if stderr_data is not None:
                    stderr_file.write(stderr_data)

        return status

    def commit(self,
               stdout_filepath=None, stderr_filepath=None,
               stdout_through=True, stderr_through=True):
        return self.commit_ver2(stdout_filepath, stderr_filepath, stderr_through, stderr_through)

    def commit_ver1(self,
                    stdout_filepath=None, stderr_filepath=None,
                    stdout_through=True, stderr_through=True):
        stdout = self._procs[-1].stdout
        stderr = self._procs[-1].stderr

        stdout_file = None
        if stdout_filepath is not None:
            stdout_file = open(stdout_filepath, mode='a')
        stderr_file = None
        if stderr_filepath is not None:
            stderr_file = open(stderr_filepath, mode='a')

        # see also.
        # http://qiita.com/FGtatsuro/items/0f68ab9c1bcad9c4b320
        # https://gist.github.com/mattbornski/3299031
        with io.open(stdout.fileno(), closefd=False) as out_stream, io.open(stderr.fileno(), closefd=False) as err_stream:
            logger.info("process commit(): stream opened.")
            for line in out_stream:
                if stdout_through:
                    sys.stdout.write(line)
                if stdout_file is not None:
                    self._nonblocking_pass(stdout_file)
                    stdout_file.write(line)
            for line in err_stream:
                if stderr_through:
                    sys.stderr.write(line)
                if stderr_file is not None:
                    stderr_file.write(line)

        self._procs[-1].wait()
        status = self._procs[-1].returncode

        for p in self._procs:
            p.stdout.close()
            p.stderr.close()
        self._procs = []

        if stdout_file is not None:
            stdout_file.close()
        if stderr_file is not None:
            stderr_file.close()

        return status

    def commit_ver2(self,
                    stdout_filepath=None, stderr_filepath=None,
                    stdout_through=True, stderr_through=True):
        """
        stderrが多すぎると刺さる。
        恐らくasyncioを使うと改善されると思われる

        Args:
            stdout_filepath ([type], optional): [description]. Defaults to None.
            stderr_filepath ([type], optional): [description]. Defaults to None.
            stdout_through (bool, optional): [description]. Defaults to True.
            stderr_through (bool, optional): [description]. Defaults to True.

        Returns:
            [type]: [description]
        """
        # logger.info("commit(): start")
        # stdout = self._procs[-1].stdout
        # stderr = self._procs[-1].stderr

        # stdout, stderrをノンブロッキングモードにする
        # self._set_nonblocking(stdout)
        # self._set_nonblocking(stderr)

        stdout_file = None
        if stdout_filepath is not None:
            stdout_file = open(stdout_filepath, mode='a')
        stderr_file = None
        if stderr_filepath is not None:
            stderr_file = open(stderr_filepath, mode='a')

        # see also.
        # http://qiita.com/FGtatsuro/items/0f68ab9c1bcad9c4b320
        # https://gist.github.com/mattbornski/3299031
        out_stream = io.open(self._procs[-1].stdout.fileno(), closefd=False)
        err_stream = None
        if self._procs[-1].stderr == subprocess.PIPE:
            err_stream = io.open(self._procs[-1].stderr.fileno(), closefd=False)

        stdout_filehandles = []
        if stdout_through:
            stdout_filehandles = [sys.stdout, stdout_file]
        else:
            stdout_filehandles = [stdout_file]

        stderr_filehandles = []
        if stderr_through:
            stderr_filehandles = [sys.stderr, stderr_file]
        else:
            stderr_filehandles = [stderr_file]

        while self._procs[-1].poll() is None:
            self._nonblocking_pass(out_stream, stdout_filehandles)
            if err_stream is not None:
                self._nonblocking_pass(err_stream, stderr_filehandles)

        self._procs[-1].wait()
        status = self._procs[-1].returncode

        self._nonblocking_pass(out_stream, stdout_filehandles, True)
        if err_stream is not None:
            self._nonblocking_pass(err_stream, stderr_filehandles, True)

        for p in self._procs:
            if p.stdout is not None:
                p.stdout.close()
            if p.stderr is not None:
                p.stderr.close()
        self._procs = []

        if stdout_file is not None:
            stdout_file.close()
        if stderr_file is not None:
            stderr_file.close()

        return status

    def _is_shell_script(self, cmd):
        """shell modeを決める
        cmdが文字列の場合はワイルドカードが使える(shell: True)。

        Args:
            cmd(str or list): [description]

        Returns:
            bool: shell modeかどうか(True)
        """        '''
        '''
        assert(isinstance(cmd, (str, list)))
        return isinstance(cmd, str)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
