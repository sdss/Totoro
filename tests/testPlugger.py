#!/usr/bin/env python
# encoding: utf-8
"""
testPlugger.py.py

Created by José Sánchez-Gallego on 14 May 2015.
Licensed under a 3-clause BSD license.

Revision history:
    14 May 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.internal.manga.Totoro.scheduler import Plugger
from sdss.internal.manga.Totoro import TotoroDBConnection
from collections import OrderedDict
import unittest

db = TotoroDBConnection()
session = db.session


class TestPlugger(unittest.TestCase):

    def setUp(self):
        """Sets default priorities for plates 8550 and 7443."""

        with session.begin(subtransactions=True):

            plate8550 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 8550).one()
            plate8550.plate_pointings[0].priority = 5

            plate7443 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 7443).one()
            plate7443.plate_pointings[0].priority = 1

    def tearDown(self):
        """Sets default priorities for plates 8550 and 7443."""

        with session.begin(subtransactions=True):

            plate8550 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 8550).one()
            plate8550.plate_pointings[0].priority = 5

            plate7443 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 7443).one()
            plate7443.plate_pointings[0].priority = 1

    def test57157(self):
        """Tests Plugger with MJD=57157."""

        plugger = Plugger(startDate=2457157.76042, endDate=2457157.95)

        validResult = OrderedDict(
            [(1, 8312), (3, 8486), (4, 8550),
             ('cart_order', [9, 8, 7, 2, 5, 6, 1, 3, 4])])

        self.assertEqual(validResult, plugger.getASOutput())

    def testForcePlug(self):
        """Tests Plugger with a plate with priotity=10."""

        # Changes the priority of plates 8550 and 7443 to 10
        with session.begin(subtransactions=True):

            plate8550 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 8550).one()
            plate8550.plate_pointings[0].priority = 10

            plate7443 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 7443).one()
            plate7443.plate_pointings[0].priority = 10

        # Now we run the plugger
        plugger = Plugger(startDate=2457185.64931, endDate=2457185.82347)

        # We expect the same result as before but with 8550 and 7443 assigned
        # first.
        validResult = OrderedDict(
            [(1, 8482), (3, 8486), (4, 7443), (5, 8550),
             ('cart_order', [9, 8, 7, 2, 6, 3, 1, 4, 5])])

        self.assertEqual(validResult, plugger.getASOutput())


if __name__ == '__main__':
    unittest.main()
