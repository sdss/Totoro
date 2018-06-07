#!/usr/bin/env python
# encoding: utf-8
"""

test_DatabaseConnect.py

Created by José Sánchez-Gallego on 10 Oct 2015.
Licensed under a 3-clause BSD license.

Revision history:
    10 Oct 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division, print_function

import configparser
import os
import textwrap
import unittest
import warnings
from builtins import str

from future import standard_library

from SDSSconnect import DatabaseConnection
from SDSSconnect.DatabaseConnection import readProfile
from SDSSconnect.exceptions import SDSSconnectError


standard_library.install_aliases()


class TestDatabaseConnect(unittest.TestCase):
    """Test suite for DatabaseConnect."""

    @classmethod
    def setUpClass(cls):
        """Sets up the test suite."""

        cls.tmpProfileSimple = os.path.join(os.path.expanduser('~'), 'test_profile_simple.ini')

        # Creates a temporary profile file.
        profileSimpleText = """
            [test]
            user: sdss
            password: sdsspass
            host: localhost
            port: 5432
            database: test

            [test2]
            user: sdssdb
            host: localhost
            port: 5432
            database: apodb
            """

        with open(cls.tmpProfileSimple, 'wb') as output:
            output.write(textwrap.dedent(profileSimpleText))

        cls.tmpProfileDefaults = os.path.join(os.path.expanduser('~'), 'test_profile_defaults.ini')

        profileDefaultsText = """
            [DEFAULT]
            user: sdssdb
            host: localhost
            port: 5432
            database: apodb

            [test]
            user: sdss
            password: sdsspass
            host: localhost
            port: 5432
            database: test
            """

        with open(cls.tmpProfileDefaults, 'wb') as output:
            output.write(textwrap.dedent(profileDefaultsText))

    @classmethod
    def tearDownClass(cls):
        """Tears down the test suite."""

        if os.path.exists(cls.tmpProfileSimple):
            os.remove(cls.tmpProfileSimple)

        if os.path.exists(cls.tmpProfileDefaults):
            os.remove(cls.tmpProfileDefaults)

    def testConfigurationFile(self):
        """Tests reading an encrypted profile."""

        cParser = configparser.ConfigParser()
        cParser.read(self.tmpProfileSimple)

        self.assertEqual(cParser.get('test', 'user'), 'sdss')
        self.assertEqual(cParser.getint('test', 'port'), 5432)

        config2 = readProfile(path=self.tmpProfileDefaults)
        self.assertEqual(config2[0]['user'], 'sdssdb')
        self.assertEqual(config2[1], 'DEFAULT')

        with warnings.catch_warnings(record=True) as ww:
            warnings.simplefilter('always')
            config3 = readProfile(path=self.tmpProfileSimple)
            self.assertIn('no default profile found. '
                          'Using first profile: test', str(ww[-1].message))

        self.assertEqual(config3[0]['user'], 'sdss')
        self.assertEqual(config3[1], 'test')

        config4 = readProfile('test2', path=self.tmpProfileSimple)
        self.assertEqual(config4[1], 'test2')
        self.assertIn('password', config4[0])
        self.assertEqual(config4[0]['password'], '')

    def testConnection(self):
        """Tests connecting to the test database."""

        self.assertEqual(len(DatabaseConnection.listConnections()), 0)

        db = DatabaseConnection('test2', profilePath=self.tmpProfileSimple)

        self.assertEqual(len(DatabaseConnection.listConnections()), 1)

        # Runs a simple query in the test DB
        session = db.Session()
        with session.begin():
            plate = session.query(db.plateDB.Plate).get(10639)

        self.assertEqual(plate.plate_id, 7495)

    def testMultipleConnections(self):
        """Tests creating multiple connections."""

        connDefault = DatabaseConnection('test2', profilePath=self.tmpProfileSimple)
        conn2 = DatabaseConnection('test2', profilePath=self.tmpProfileSimple)
        self.assertIs(connDefault, conn2)
        self.assertEqual(connDefault.profile, 'test2')

        with warnings.catch_warnings(record=True) as ww:
            warnings.simplefilter('always')
            conn3 = DatabaseConnection('test2', new=True, profilePath=self.tmpProfileSimple)
            self.assertIn('overwritting profile test2', str(ww[-1].message))

        self.assertIs(connDefault, conn2)
        self.assertIsNot(connDefault, conn3)
        self.assertIsNot(conn2, conn3)

        conn4 = DatabaseConnection(
            'test2', new=True, name='testConn', default=True, profilePath=self.tmpProfileSimple)
        self.assertItemsEqual(conn4.listConnections(), ['test2', 'testConn'])
        self.assertItemsEqual(conn4.listConnections(), DatabaseConnection.listConnections())
        self.assertIsNot(conn4, conn3)
        self.assertIsNot(conn4, connDefault)
        self.assertEqual(conn4.getDefaultConnectionName(), 'testConn')

        conn5 = DatabaseConnection(profilePath=self.tmpProfileDefaults)
        self.assertIs(conn5, conn4)

        with self.assertRaises(SDSSconnectError):
            DatabaseConnection('production')


if __name__ == '__main__':
    unittest.main()
