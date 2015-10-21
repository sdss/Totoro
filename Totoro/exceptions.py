# Exceptions for pyPhotom classes and functions

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


class TotoroError(Exception):
    """Base exception for Totoro. Other exceptions should inherit this."""
    pass


class TotoroWarning(Warning):
    """Base warning for Totoro."""


class TotoroUserWarning(TotoroWarning):
    """The primary warning class."""
