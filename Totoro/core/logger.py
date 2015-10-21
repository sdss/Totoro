# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines a logging class based on the built-in logging module.
The module is heavily based on the astropy logging system.
"""

from __future__ import print_function

import os
import sys
import logging
from logging import FileHandler
import warnings
from . import colourPrint
import shutil
from .defaults import *

# Initialize by calling initLog()
log = None


def initLog():

    logging.setLoggerClass(TotoroLogger)
    log = logging.getLogger('Totoro')
    log._set_defaults()

    return log


class MyFormatter(logging.Formatter):

    warning_fmp = '%(asctime)s - %(levelname)s: %(message)s [%(origin)s]'
    info_fmt = '%(asctime)s - %(levelname)s: %(message)s'

    def __init__(self, fmt='%(levelno)s: %(msg)s'):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._fmt = MyFormatter.info_fmt

        elif record.levelno == logging.INFO:
            self._fmt = MyFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._fmt = MyFormatter.info_fmt

        elif record.levelno == logging.WARNING:
            self._fmt = MyFormatter.warning_fmp

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        return result


Logger = logging.getLoggerClass()
fmt = MyFormatter()


class TotoroLogger(Logger):
    """This class is used to set up the logging system.

    The main functionality added by this class over the built-in
    logging.Logger class is the ability to keep track of the origin of the
    messages, the ability to enable logging of warnings.warn calls and
    exceptions, and the addition of colorized output and context managers to
    easily capture messages to a file or list.

    """

    def saveLog(self, path):
        shutil.copyfile(self.logFilename, os.path.expanduser(path))

    def _showwarning(self, *args, **kwargs):

        warning = args[0]
        message = '{0}: {1}'.format(warning.__class__.__name__, args[0])
        mod_path = args[2]

        mod_name = None
        mod_path, ext = os.path.splitext(mod_path)
        for name, mod in sys.modules.items():
            path = os.path.splitext(getattr(mod, '__file__', ''))[0]
            if path == mod_path:
                mod_name = mod.__name__
                break

        if mod_name is not None:
            self.warning(message, extra={'origin': mod_name})
        else:
            self.warning(message)

    def _stream_formatter(self, record):
        """The formatter for standard output."""
        if record.levelno < logging.DEBUG:
            print(record.levelname, end='')
        elif(record.levelno < logging.INFO):
            colourPrint(record.levelname, 'magenta', end='')
        elif(record.levelno < logging.WARN):
            colourPrint(record.levelname, 'green', end='')
        elif(record.levelno < logging.ERROR):
            colourPrint(record.levelname, 'brown', end='')
        else:
            colourPrint(record.levelname, 'red', end='')

        if not hasattr(record, 'origin') or record.origin == '':
            record.message = '{0}'.format(record.msg)
        else:
            record.message = '{0} [{1:s}]'.format(record.msg, record.origin)

        print(': ' + record.message)

    def _set_defaults(self):
        """Reset logger to its initial state."""

        # Remove all previous handlers
        for handler in self.handlers[:]:
            self.removeHandler(handler)

        # Set levels
        self.setLevel('DEBUG')

        # Set up the stdout handler
        sh = logging.StreamHandler()
        sh.emit = self._stream_formatter
        self.addHandler(sh)
        sh.setLevel(LOG_LEVEL())

        # Set up the main log file handler if requested (but this might fail if
        # configuration directory or log file is not writeable).

        log_file_path = LOG_FILE_PATH()
        if not os.path.exists(LOG_DIR()):
            os.mkdir(LOG_DIR())

        try:
            if log_file_path == '':
                log_file_path = os.path.join(
                    LOG_DIR(), 'totoro.log')
            else:
                log_file_path = os.path.expanduser(log_file_path)

            fh = FileHandler(log_file_path, mode=LOG_FILE_MODE())
        except (IOError, OSError) as e:
            warnings.warn(
                'log file {0!r} could not be opened for writing: '
                '{1}'.format(log_file_path, unicode(e)), RuntimeWarning)
        else:
            fh.setFormatter(fmt)
            fh.setLevel(LOG_FILE_LEVEL())
            self.addHandler(fh)

        self.logFilename = log_file_path
        warnings.showwarning = self._showwarning
