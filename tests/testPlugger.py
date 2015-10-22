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
from sdss.internal.manga.Totoro import TotoroDBConnection, config
from collections import OrderedDict
import unittest

db = TotoroDBConnection()
session = db.session


class TestPlugger(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Sets up the test suite."""

        cls.originalOfflineCarts = config['offlineCarts']

        config['offlineCarts'] = [2]

        with session.begin():

            # Restores priorities
            plate8550 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 8550).one()
            plate8550.plate_pointings[0].priority = 5

            plate7443 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 7443).one()
            plate7443.plate_pointings[0].priority = 1

    def tearDown(self):
        """Restores changes made to the DB during the tests."""

        with session.begin():

            # Restores priorities
            plate8550 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 8550).one()
            plate8550.plate_pointings[0].priority = 5

            plate7443 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 7443).one()
            plate7443.plate_pointings[0].priority = 1

            # Restores the list of active pluggings
            for ii in range(1, 4):
                session.delete(
                    session.query(db.plateDB.ActivePlugging).get(ii))

            session.add(db.plateDB.ActivePlugging(plugging_pk=70172, pk=1))
            session.add(db.plateDB.ActivePlugging(plugging_pk=70173, pk=2))
            session.add(db.plateDB.ActivePlugging(plugging_pk=70124, pk=3))

        # Restores offline carts
        config['offlineCarts'] = self.originalOfflineCarts

    def test57157(self):
        """Tests Plugger with MJD=57157."""

        plugger = Plugger(startDate=2457157.76042, endDate=2457157.95)

        validResult = OrderedDict(
            [(1, 8312), (3, 8486), (4, 8550),
             ('cart_order', [9, 8, 7, 2, 5, 6, 3, 1, 4])])

        self.assertEqual(validResult, plugger.getASOutput())

    def testForcePlug(self):
        """Tests Plugger with a plate with priotity=10."""

        # Changes the priority of plates 8550 and 7443 to 10
        with session.begin():

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

    def testCartOrder(self):
        """Tests the cart order returned by the Plugger."""

        # Sets a custom list of active pluggings.
        with session.begin():

            for ii in range(1, 4):
                session.delete(
                    session.query(db.plateDB.ActivePlugging).get(ii))

            session.add(db.plateDB.ActivePlugging(plugging_pk=70003, pk=1))
            session.add(db.plateDB.ActivePlugging(plugging_pk=68904, pk=2))
            session.add(db.plateDB.ActivePlugging(plugging_pk=70124, pk=3))

        # Calls the Plugger without dates
        plugger = Plugger(startDate=None, endDate=None)

        self.assertEqual(plugger.getASOutput()['cart_order'],
                         [9, 8, 7, 4, 5, 6, 3, 2, 1])

    def testOfflineCarts(self):
        """Tests if Plugger works if a plate is plugged in an offline cart."""

        # Makes cart 3 offline
        config['offlineCarts'] = [3]

        # Runs Plugger
        plugger = Plugger(startDate=2457182.64792, endDate=2457182.79097)
        output = plugger.getASOutput()

        self.assertEqual(output,
                         OrderedDict([(1, 8482), (2, 8312), (3, 8486),
                                      ('cart_order', [9, 8, 7, 4, 5, 6,
                                                      1, 3, 2])]))

if __name__ == '__main__':
    unittest.main()
