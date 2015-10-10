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

from __future__ import division
from __future__ import print_function
import unittest
import configparser
import os
from simplecrypt import encrypt, decrypt
import StringIO
from SDSSconnect import DatabaseConnection


class TestDatabaseConnect(unittest.TestCase):
    """Test suite for DatabaseConnect."""

    @classmethod
    def setUpClass(cls):
        """Sets up the test suite."""

        cls.password = 'secure'

        cls.tmpProfiles = os.path.join(os.path.expanduser('~'),
                                       'test_profiles.ini')

        # Creates a temporary profile file.
        profileText = """
            [test]
            user: sdss
            password: sdsspass
            host: localhost
            port: 5432
            database: test
            """

        with open(cls.tmpProfiles, 'wb') as output:
            ciphertext = encrypt(cls.password, profileText)
            output.write(ciphertext)

    @classmethod
    def tearDownClass(cls):
        """Tears down the test suite."""

        if os.path.exists(cls.tmpProfiles):
            os.remove(cls.tmpProfiles)

    def testConfigurationFile(self):
        """Tests reading an encrypted profile."""

        plainText = decrypt(self.password,
                            open(self.tmpProfiles, 'r').read()).decode('utf8')
        buf = StringIO.StringIO(plainText)

        cParser = configparser.ConfigParser()
        cParser.readfp(buf)

        self.assertEqual(cParser.get('test', 'user'), 'sdss')
        self.assertEqual(cParser.getint('test', 'port'), 5432)

    def testConnection(self):
        """Tests connecting to the test database."""

        self.assertEqual(len(DatabaseConnection.listConnections()), 0)

        db = DatabaseConnection('test')

        self.assertEqual(len(DatabaseConnection.listConnections()), 1)

        # Runs a simple query in the test DB
        session = db.Session()
        with session.begin():
            actPlug = session.query(db.plateDB.ActivePlugging).get(1)

        self.assertEqual(actPlug.plugging_pk, 70172)

    def testMultipleConnections(self):
        """Tests creating multiple connections."""

        connDefault = DatabaseConnection('test')
        conn2 = DatabaseConnection('test')
        self.assertIs(connDefault, conn2)

        conn3 = DatabaseConnection('test', new=True)
        self.assertIs(connDefault, conn2)
        self.assertIsNot(connDefault, conn3)
        self.assertIsNot(conn2, conn3)

        conn4 = DatabaseConnection('test', new=True, name='testConn')
        self.assertItemsEqual(conn4.listConnections(), ['default', 'testConn'])
        self.assertItemsEqual(conn4.listConnections(),
                              DatabaseConnection.listConnections())
        self.assertIsNot(conn4, conn3)
        self.assertIsNot(conn4, connDefault)


if __name__ == '__main__':
    unittest.main()
