#!/usr/bin/env python
# encoding: utf-8
"""
apoDB.py

Created by José Sánchez-Gallego on 30 Mar 2014.
Licensed under a 3-clause BSD license.

Sets up the APO database configuration and creates the objects necessary to
connect to plateDB and mangaDB.

Revision history:
    30 Mar 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import sys
from sdss.internal.manga.Totoro import config, log

try:
    from sdss.internal.database.DatabaseConnection import DatabaseConnection
except ImportError as e:
    print('Couldn\'t find DatabaseConnection:', e)
    sys.exit(1)


class Singleton(type):

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(
                Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class TotoroDBConnection(object):

    __metaclass__ = Singleton

    def __init__(self):

        # The database connection string depends on the version of SQLAlchemy.
        self.databaseConnectionString = (
            'postgresql+psycopg2://{user}:{password}@{host}:{port}/{database}'
            .format(**config['dbConnection']))

        log.debug('Creating connection with DB.')

        # Intial database connection creation and instances to be exported.
        self.db = DatabaseConnection(
            database_connection_string=self.databaseConnectionString,
            expire_on_commit=False)

        self.engine = self.db.engine
        self.metadata = self.db.metadata
        self.Session = self.db.Session
        self.Base = self.db.Base

        import sdss.internal.database.apo.platedb.ModelClasses as plateDB
        import sdss.internal.database.apo.mangadb.ModelClasses as mangaDB

        self.plateDB = plateDB
        self.mangaDB = mangaDB
        self.session = self.Session()
