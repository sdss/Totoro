#!/usr/bin/env python
# encoding: utf-8
"""

exceptions.py

Created by José Sánchez-Gallego on 15 Oct 2015.
Licensed under a 3-clause BSD license.

Revision history:
    15 Oct 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division, print_function, unicode_literals


class SDSSconnectError(Exception):
    """Base exception for SDSSconnect. Other exceptions should inherit this."""
    pass


class SDSSconnectWarning(Warning):
    """Base warning for SDSSconnect."""
    pass


class SDSSconnectUserWarning(UserWarning, SDSSconnectWarning):
    """The primary warning class."""
    pass
