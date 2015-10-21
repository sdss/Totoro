#!/usr/bin/env python
# encoding: utf-8
"""
connection.py

Created by José Sánchez-Gallego on 30 Oct 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

import os
from ..core import ConfigObject
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base

# APO database configuration parameters
_db_config = {
    'user': 'albireo',
    'password': '',
    'database': 'apo_platedb',
    'host': 'localhost',
    'port': 5432
}

# These value can be overridden with environment variables
if 'PLATEDB_USER' in os.environ:
    _db_config['user'] = os.environ['PLATEDB_USER']
if 'PLATEDB_PASSWORD' in os.environ:
    _db_config['password'] = os.environ['PLATEDB_PASSWORD']
if 'PLATEDB_HOST' in os.environ:
    _db_config['host'] = os.environ['PLATEDB_HOST']
if 'PLATEDB_PORT' in os.environ:
    _db_config['port'] = os.environ['PLATEDB_PORT']
if 'PLATEDB_DB' in os.environ:
    _db_config['database'] = os.environ['PLATEDB_DB']

# Finally, the options can be overridden if a configuration file
# exists.

USER = ConfigObject('platedb_user', _db_config['user'],
                    'The user for plateDB', 'plateDB')
PASSWORD = ConfigObject('platedb_password', _db_config['password'],
                        'The password for plateDB', 'plateDB')
HOST = ConfigObject('platedb_host', _db_config['host'],
                    'The host for plateDB', 'plateDB')
PORT = ConfigObject('platedb_port', _db_config['port'],
                    'The port for plateDB', 'plateDB')
DB_NAME = ConfigObject('platedb_db_name', _db_config['database'],
                       'The database name for plateDB', 'plateDB')


class Singleton(type):
    def __init__(cls, name, bases, dict):
        super(Singleton, cls).__init__(name, bases, dict)
        cls.instance = None

    def __call__(cls, *args, **kw):
        if cls.instance is None:
            cls.instance = super(Singleton, cls).__call__(*args, **kw)
        return cls.instance


class DatabaseConnection(object):
    """Creates a connection to a database.

    This class is a singleton to avoid that different instances
    access the database simultaneously.
    """

    __metaclass__ = Singleton

    def __init__(self, connectionString=None):

        if connectionString is None:
            connectionString = 'postgresql+psycopg2://' + \
                '{0}:{1}@{2}:{3}/{4}'.format(
                    USER(), PASSWORD(), HOST(), PORT(), DB_NAME())

        self.engine = create_engine(connectionString, echo=False)
        self.metadata = MetaData()
        Session = sessionmaker(bind=self.engine)
        self.Base = declarative_base(bind=self.engine)
        self.session = Session()

    def close(self):
        self.session.close()
