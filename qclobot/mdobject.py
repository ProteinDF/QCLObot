#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .taskobject import TaskObject

class MdObject(TaskObject):
    """Operate MD task

    : Example
    >>> md1 = MdObject(name='md_test1')
    """

    def __init__(self, *args, **kwargs):
        super(MdObject, self).__init__(*args, **kwargs)

    def opt(self):
        pass

if __name__ == '__main__':
    import doctest
    doctest.testmod()
