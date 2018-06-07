#!/usr/bin/env python
# encoding: utf-8
"""

db.py

Created by José Sánchez-Gallego on 25 Oct 2015.
Licensed under a 3-clause BSD license.

Revision history:
    25 Oct 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division, print_function

from SDSSconnect import DatabaseConnection
from Totoro import config


def getConfigurationProfiles():
    """Returns a dictionary with all currently configured DB profiles."""

    profiles = {}
    for kk in config:
        if 'dbConnection' in kk and kk != 'dbConnection':
            profileName = config[kk]['name'].lower()
            profiles[profileName] = config[kk]
            if 'password' not in profiles[profileName]:
                profiles[profileName]['password'] = ''

    return profiles


def getConnection(profile=None):
    """Returns a connection.

    If `profile=None`, the default connection is returned."""

    # To avoid circular import errors
    from Totoro.utils.utils import checkOpenSession

    configProfiles = getConfigurationProfiles()

    if len(DatabaseConnection.listConnections()) > 0 and profile is None:
        return DatabaseConnection.getDefaultConnection()
    elif len(DatabaseConnection.listConnections()) == 0 and profile is None:
        # Creates the default DB connection
        databaseConnectionString = (
            'postgresql+psycopg2://{user}:{password}@{host}:{port}/{database}'
            .format(**config['dbConnection']))
        dbConn = DatabaseConnection(
            databaseConnectionString=databaseConnectionString,
            new=True,
            name=config['dbConnection']['name'],
            default=True)
        checkOpenSession()
        return dbConn
    else:
        if profile.lower() in DatabaseConnection.listConnections():
            return DatabaseConnection.getConnection(profile.lower())
        else:
            if profile.lower() in configProfiles:
                databaseConnectionString = ('postgresql+psycopg2://{user}:{password}@'
                                            '{host}:{port}/{database}'
                                            .format(**configProfiles[profile.lower()]))
                dbConn = DatabaseConnection(
                    databaseConnectionString=databaseConnectionString,
                    new=True,
                    name=profile.lower())
                checkOpenSession()
                return dbConn
            else:
                raise ValueError('profile {0} does not exist'.format(profile))


def getConnectionFull(profile=None):
    """Returns a connection, its session, plateDB and mangaDB."""

    dbConn = getConnection(profile=profile)
    return dbConn, dbConn.Session, dbConn.plateDB, dbConn.mangaDB


def setDefaulProfile(profile):
    """Sets a profile as default."""

    if len(DatabaseConnection.listConnections()) > 0:
        if DatabaseConnection.getDefaultConnectionName() == 'profile':
            return

    if profile not in getConfigurationProfiles():
        raise ValueError('profile {0} does not exist'.format(profile))

    if profile in DatabaseConnection.listConnections():
        DatabaseConnection.setDefaultConnection(profile)
    else:
        db = getConnection(profile=profile)
        db.setDefaultConnection(profile)
