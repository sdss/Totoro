#!/usr/bin/env python
# encoding: utf-8
"""
configuration.py

Created by José Sánchez-Gallego on 10 Dec 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

Includes classes and functions to read and write the configuration
file for Totoro. Uses a simplified version of the astropy configuration
system.

"""

import ConfigParser
import os
from types import NoneType


TOTORO_CONFIG_PATH = os.path.join(os.getenv('HOME'), '.totoro/totoro.config')


class ConfigObject(object):
    """Similar to ConfigurationItem() in astropy but simpler.

    Defines a configuration object with a default value. The value
    is returned by calling the object. If the configuration item exists
    in the configuration file, returns that value. Otherwise, returns
    the default value.

    Parameters
    ----------

    name : str
        The name of the configuration option.
    defaultValue
        The value to use if no configuration file is present.
    comment : str, optional
        A comment to be saved with the option.
    section : str, optional
        The section of the configuration file to save the file. The
        default is Totoro.

    """

    def __init__(self, name, defaultValue, comment='', section='Totoro'):

        self.defaultValue = self._evalValue(defaultValue)
        self.name = name
        self.comment = comment
        self.section = section
        self.value = None

        self.configPath = os.path.dirname(TOTORO_CONFIG_PATH)

        self._configparser = self._readConfigParser()

    def __call__(self):
        """Returns the value in the configuration file or the default."""

        if self.value is not None:
            return self.value

        try:
            return self._evalValue(
                self._configparser.get(
                    self.section, self.name))
        except:
            return self.defaultValue

    def save(self):
        """Saves the default value to the configuration file."""

        if self._configparser is None:
            self._createConfigParser()

        if not self._configparser.has_section(self.section):
            self._configparser.add_section(self.section)

        if self.comment == '':
            self._configparser.set(self.section, self.defaultValue)
        else:
            self._configparser.set(self.section, self.defaultValue +
                                   '   ; {0}'.format(self.comment))

    def _createConfigParser(self):
        self._configparser = ConfigParser.SafeConfigParser()
        self._configparser.optionxform = str

    def _evalValue(self, value):
        """Evaluates the type of a value."""

        if isinstance(value, NoneType):
            return None
        try:
            return eval(value)
        except:
            return value

    def _readConfigParser(self):
        """Reads Totoro configuration file, if any."""

        if os.path.exists('./totoro.config'):
            configFile = './totoro.config'
        elif os.path.exists(TOTORO_CONFIG_PATH):
            configFile = TOTORO_CONFIG_PATH
        else:
            return None

        config = ConfigParser.SafeConfigParser()
        config.optionxform = str
        config.read(configFile)

        return config
