#!/usr/bin/env python
# -*- coding: utf-8 -*-

import qclobot as qclo


def main():
    p = qclo.process.Process()

    p.cmd("ls -al")
    p.pipe("grep pdb")

    code = p.commit2("stdout.txt", "stderr.txt")
    print(code)


if __name__ == '__main__':
    main()
