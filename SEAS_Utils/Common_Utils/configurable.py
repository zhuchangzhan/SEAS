#!/usr/bin/env python
#
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Base class for configurable objects


This code is inherited the configurable.py in TSIG code written by
Martin Owens, et al. 
"""

import inspect
import configobj
import os

class Configuration(configobj.ConfigObj):
    def __init__(self, filename=None):
        """
        Create configuration, from file if specified.
        """
        if filename:
            if not os.path.exists(filename):
                raise Exception("No configuration found at %s" % filename)
        super(Configuration, self).__init__(filename)


class MetaObject(type):
    def __call__(self, cls, *args, **kwargs):
        if args and isinstance(args[0], dict):
            (config, args) = (args[0], args[1:])
            # Combine all the configurations together into one kwargs
            for key in cls.config_keys or [cls.__name__]:
                kwargs.update(config.get(key, {}))

            # Protect classes that won't want to accept **kwargs by removing
            # Any keys that don't match known arguments to the class.
            spec = inspect.getargspec(cls.__init__)
            if not spec.keywords:
                for key in list(kwargs):
                    if key not in spec.args:
                        kwargs.pop(key)

        obj = cls.__new__(cls, *args, **kwargs)
        obj.__init__(*args, **kwargs)
        return obj


class ConfigurableObject(object):
    """
    When constructing an child class, you can specify a list of config_keys
    which are strings used on the first level of the configuration structure.

    For example:

       class A(ConfigurableObject):
           config_keys = ['B', 'C']

    Class 'A' in this case would capture arguments from

       {'A': {...}, 'B': {...'}}

    in that order. So any duplicates would clober the previous one as the new
    keyword argument dictionary is built.

    The default is the classes own name, so 'A' in this case.
    """
    __metaclass__ = MetaObject
    config_keys = None


if __name__ == '__main__':
    # Some basic tests
    class B(ConfigurableObject):
        config_keys = ['C', 'A']
        def __init__(self, a=4, b=10):
            self.a = a
            self.b = b

    assert B().a == 4
    assert B(a=10).a == 10
    assert B({}).a == 4
    assert B({'C': {'a': 20, 'd': False}}).a == 20
    assert B({'C': {'a': 20}, 'A': {'a': 9}}).a == 9
    print "OK"
