#!/usr/bin/env python
# encoding: utf-8
"""
readPath.py

Created by José Sánchez-Gallego on 25 Apr 2014.
Licensed under a 3-clause BSD license.

Revision history:
    25 Apr 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import os


def readPath(path):
    """Get a Totoro-formatted path and returns the real path.

    Paths are expanded depending on the first characher of the input string.
    If the first character is '+', the path is considered to be
    Totoro-internal relative to the root of the package. Otherwise, the path
    is expanded using os.path.expandvars and os.path.expanduser. So,
    environment variables and user tildes '~' are valid in any path.

    """

    if path[0] == '+':
        return os.path.realpath(
            os.path.join(
                os.path.dirname(__file__), path[1:]))

    else:
        return os.path.realpath(os.path.expanduser(os.path.expandvars(path)))
