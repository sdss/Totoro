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

import yaml
import os
from Totoro import __TOTORO_CONFIG_PATH__, __DEFAULT_CONFIG_FILE__
from Totoro.exceptions import TotoroError


def getConfiguration(default=None):

    from .. import __DEFAULT_CONFIG_FILE__, __TOTORO_CONFIG_PATH__

    if default is None:
        default = __DEFAULT_CONFIG_FILE__

    config = TotoroConfig(default)

    if os.path.exists(__TOTORO_CONFIG_PATH__):
        # If a custom configuration file exists, updates default values.
        config.updateFromFile(__TOTORO_CONFIG_PATH__)

    return config


class TotoroConfig(dict):

    def __init__(self, configurationFile):

        if not os.path.exists(configurationFile):
            raise TotoroError('configuration file', configurationFile,
                              'not found.')

        self._rawData = open(configurationFile).read()
        self._initFromRaw()

    def _initFromRaw(self):

        yamlData = yaml.load(self._rawData)

        if yamlData is None:
            yamlData = {}

        dict.__init__(self, yamlData)

        self._checkDBConnection()

    def save(self, path=__TOTORO_CONFIG_PATH__):
        outUnit = open(path, 'w')
        yaml.dump(dict(self), outUnit, default_flow_style=False)
        outUnit.close()

    def createTemplate(self, path=__TOTORO_CONFIG_PATH__):

        defaultConfig = TotoroConfig(__DEFAULT_CONFIG_FILE__)
        defaultConfig.save()

    def updateFromFile(self, file):

        newRawData = open(file).read()
        self._rawData += newRawData
        self._initFromRaw()

    def _checkDBConnection(self):

        if 'TOTORO_DB_CONNECTION' in os.environ:
            result = self._assignDBConnection(
                os.environ['TOTORO_DB_CONNECTION'])
            if result is True:
                return

        if isinstance(self['dbConnection'], dict):
            pass
        elif isinstance(self['dbConnection'], basestring):
            result = self._assignDBConnection(self['dbConnection'])
            if result is False:
                raise TotoroError('dbConnection could not be configured. '
                                  'Please, check your totoro.yaml and '
                                  '$TOTORO_DB_CONNECTION')

    def _assignDBConnection(self, value):

        if value.lower() in ['production', 'tunnel', 'dev', 'local']:
            self['dbConnection'] = self['dbConnection' + value.title()]
            return True
        else:
            return False
