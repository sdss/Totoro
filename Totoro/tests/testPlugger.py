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
from Totoro.scheduler import Plugger
from Totoro import config
from Totoro.db import getConnection
from collections import OrderedDict
import unittest

db = getConnection('test')
session = db.Session()


class TestPlugger(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Sets up the test suite."""

        cls.originalOfflineCarts = config['offlineCarts']
        cls._restoreActivePluggings()

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

            plate8081 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 8081).one()
            plate8081.plate_location_pk = 28

        self._restoreActivePluggings()

        # Restores offline carts
        config['offlineCarts'] = self.originalOfflineCarts

    @staticmethod
    def _restoreActivePluggings():
        """Restores the active pluggings for the test."""

        # Restores the list of active pluggings
        for ii in range(1, 4):
            aP = session.query(db.plateDB.ActivePlugging).get(ii)
            if aP is not None:
                session.delete(aP)

        session.add(db.plateDB.ActivePlugging(plugging_pk=70172, pk=1))
        session.add(db.plateDB.ActivePlugging(plugging_pk=70173, pk=2))
        session.add(db.plateDB.ActivePlugging(plugging_pk=70124, pk=3))

    def test57157(self):
        """Tests Plugger with MJD=57157."""

        plugger = Plugger(startDate=2457157.76042, endDate=2457157.95)

        validResult = OrderedDict(
            [(1, 8312), (3, 8486), (4, 8550),
             ('cart_order', [9, 8, 7, 5, 6, 2, 3, 1, 4])])

        self.assertEqual(validResult, plugger.getASOutput())

    def test57307(self):
        """Tests replugging when cart is offline."""

        # Moves plate 8081 to APO
        with session.begin():
            plate8081 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 8081).one()
            plate8081.plate_location_pk = 23

        plugger = Plugger(startDate=2457307.806736, endDate=2457307.998611)

        validResult = OrderedDict([(1, 8570), (3, 8486), (4, 8081), (5, 8566),
                                   ('cart_order', [9, 8, 7, 6, 2, 3, 4, 5, 1])
                                   ])

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
            [(1, 8482), (3, 8486), (4, 7443), (5, 8550), (6, 8312),
             ('cart_order', [9, 8, 7, 2, 3, 1, 6, 4, 5])])

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

        # Now we want to test if cart 2 (unplugged) is returned with high
        # cart order (but still fewer than cart 1, which has incomplete sets)
        # if we are running in a non-MaNGA night (the remaining tests already
        # check this for MaNGA nights).

        # We first unplug cart 2 and make it offline
        with session.begin():
            session.delete(session.query(db.plateDB.ActivePlugging).get(2))

        config['offlineCarts'] = [2]

        # We run the plugger for a non-MaNGA night
        pluggerNoMaNGA = Plugger(startDate=None, endDate=None)
        self.assertEqual(pluggerNoMaNGA.getASOutput()['cart_order'],
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
