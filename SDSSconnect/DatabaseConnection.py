#!/usr/bin/env python
# encoding: utf-8
"""

DatabaseConnection.py

Created by José Sánchez-Gallego on 10 Oct 2015.
Licensed under a 3-clause BSD license.

Revision history:
    10 Oct 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sqlalchemy import create_engine, MetaData
import sqlalchemy
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.event import listen
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.pool import Pool
from SDSSconnect.models import createRelationships
from SDSSconnect.models.utils import ModelWrapper, cameliseClassname
from SDSSconnect.models.utils import nullifyRelationship
from SDSSconnect.exceptions import SDSSconnectUserWarning, SDSSconnectError

import warnings
import ConfigParser
import os


__MODELS__ = ['plateDB', 'mangaDB']


def readProfile(profile=None, path=None):
    """Reads a profile and creates the appropriate connection string."""

    profilesPath = os.path.join(os.path.expanduser('~'), '.sdssconnect',
                                'profiles.ini') if path is None else path

    if not os.path.exists(profilesPath):
        raise RuntimeError('profile not found in {0}'.format(profilesPath))

    config = ConfigParser.ConfigParser()
    config.read(profilesPath)

    # If no profile is defined, we try to return the DEFAULTS section and, if
    # that fails, the first section in the profiles file.
    if profile is None:
        if len(config.defaults()) > 0:
            return (config.defaults(), 'DEFAULT')
        else:
            if len(config.sections()) == 0:
                raise RuntimeError('no sections found in {0}'
                                   .format(profilesPath))
            section = config.sections()[0]
            warnings.warn('no default profile found. Using first profile: {}'
                          .format(section), SDSSconnectUserWarning)
            returnDict, returnProfile = dict(config.items(section)), section
    else:
        if not config.has_section(profile.lower()):
            raise ValueError('profile {0} does not exist'
                             .format(profile.lower()))
        returnDict = dict(config.items(profile.lower()))
        returnProfile = profile.lower()

    # If the password is not present, we assume is empty (this is useful if we
    # don't define a password because it's set in pgpass)
    if 'password' not in returnDict:
        returnDict['password'] = ''

    return returnDict, returnProfile


def clearSearchPathCallback(dbapi_con, connection_record):
    """
    When creating relationships across schema, SQLAlchemy
    has problems when you explicitly declare the schema in
    ModelClasses and it is found in search_path.

    The solution is to set the search_path to "$user" for
    the life of any connection to the database. Since there
    is no (or shouldn't be!) schema with the same name
    as the user, this effectively makes it blank.

    This callback function is called for every database connection.

    For the full details of this issue, see:
    http://groups.google.com/group/sqlalchemy/browse_thread/
    thread/88b5cc5c12246220

    dbapi_con - type: psycopg2._psycopg.connection
    connection_record - type: sqlalchemy.pool._ConnectionRecord
    """

    cursor = dbapi_con.cursor()
    cursor.execute('SET search_path TO "$user",functions')
    dbapi_con.commit()


listen(Pool, 'connect', clearSearchPathCallback)


class DatabaseConnection(object):

    _singletons = dict()
    _defaultConnectionProfile = None

    def __new__(cls, *args, **kwargs):

        if len(args) == 1:
            kwargs['profile'] = args[0]
        elif len(args) > 1:
            raise ValueError('SDSSDatabaseConnection only '
                             'accepts one argument')

        new = kwargs.get('new', False)
        default = kwargs.get('default', False)
        profile = kwargs.get('profile', None)

        if len(cls._singletons.keys()) == 0 or new:

            newInstance = cls._createNewInstance(**kwargs)
            profile = newInstance.profile

            if profile in cls.listConnections():
                warnings.warn('overwritting profile {0}'.format(profile),
                              SDSSconnectUserWarning)

            if default or len(cls._singletons.keys()) == 0:
                cls._defaultConnectionProfile = newInstance.profile

            cls._singletons[newInstance.profile] = newInstance

            return cls._singletons[newInstance.profile]

        else:

            if profile is None and len(cls.listConnections()) == 0:
                raise SDSSconnectError('no available connection exists')
            else:
                if profile is None:
                    return cls.getDefaultConnection()
                elif profile not in cls.listConnections():
                    raise SDSSconnectError('profile {0} not in the list of '
                                           'connections. Use new=True if you '
                                           'want to create a new connection.'
                                           .format(profile))
                else:
                    return cls._singletons[profile]

    @classmethod
    def _createNewInstance(cls, profile=None, databaseConnectionString=None,
                           expireOnCommit=True, models='all',
                           profilePath=None, name=None, **kwargs):
        """Creates a new instance of the connection."""

        me = object.__new__(cls)

        if databaseConnectionString is not None:
            me.databaseConnectionString = databaseConnectionString
            me.profile = 'connection_string' if name is None else name
        else:
            profileDict, profileName = readProfile(profile=profile,
                                                   path=profilePath)
            me.databaseConnectionString = (
                'postgresql+psycopg2://{user}:{password}@'
                '{host}:{port}/{database}'.format(**profileDict))

            me.profile = profileName if name is None else name

        me.engine = create_engine(me.databaseConnectionString, echo=False)

        me.metadata = MetaData()
        me.Session = scoped_session(
            sessionmaker(bind=me.engine, autocommit=True,
                         expire_on_commit=expireOnCommit))

        me.addModels(models)

        return me

    @classmethod
    def getConnection(cls, connectionProfile):
        """Returns a named connection."""

        if connectionProfile not in cls._singletons:
            raise KeyError('connection named {0} does not exist'
                           .format(connectionProfile))

        return cls._singletons[connectionProfile]

    @classmethod
    def listConnections(cls):
        """Returns a list of all available connections."""

        return cls._singletons.keys()

    @classmethod
    def getDefaultConnectionName(cls):
        """Returns the default connection name."""

        if len(cls._singletons) == 0:
            raise SDSSconnectError('no connections available')
        elif cls._defaultConnectionProfile not in cls._singletons:
            raise SDSSconnectError('default profile {0} is not in the list of '
                                   'conections'.format(cls._defaultConnection))

        return cls._defaultConnectionProfile

    @classmethod
    def getDefaultConnection(cls):
        """Returns the default connection."""

        return cls._singletons[cls.getDefaultConnectionName()]

    @classmethod
    def setDefaultConnection(cls, connectionProfile):
        """Sets the default connection."""

        if connectionProfile not in cls._singletons:
            raise KeyError('connection profile {0} does not exist'
                           .format(connectionProfile))

        cls._defaultConnectionProfile = connectionProfile

    def addModels(self, models, overwrite=False):
        """Adds a list of model classes to the object."""

        assert isinstance(models, (list, tuple, basestring)), \
            ('models must be \'all\' or a list of strings.')

        self.metadata = MetaData(bind=self.engine)

        if models == 'all':
            models = __MODELS__

        for model in models:
            if hasattr(self, model) and not overwrite:
                continue

            self.metadata.reflect(schema=model.lower())

        self.Base = automap_base(metadata=self.metadata)
        self.Base.prepare(engine=self.engine,
                          classname_for_table=cameliseClassname)

        # We are going to create our own relationships, so we remove these.
        for model in self.Base.classes:
            for relationship in model.__mapper__.relationships.keys():
                delattr(model, str(relationship))

        createRelationships(self.Base)
        sqlalchemy.orm.configure_mappers()

        # Now that all classes have been added to the Base, we wrap them
        # depending on their schema.
        for model in models:
            if hasattr(self, model) and not overwrite:
                continue
            if model.lower() == 'platedb':
                self.plateDB = ModelWrapper(self.Base, 'Platedb_')
            elif model.lower() == 'mangadb':
                self.mangaDB = ModelWrapper(self.Base, 'Mangadb_')
