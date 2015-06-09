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

    def setUp(self):
        """Sets default priorities for plates 8550 and 7443."""

        self.originalOfflineCarts = config['offlineCarts']

        with session.begin(subtransactions=True):

            # Saves the list of active pluggings
            self.originalActivePluggings = session.query(
                db.plateDB.ActivePlugging).all()

            plate8550 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 8550).one()
            plate8550.plate_pointings[0].priority = 5

            plate7443 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 7443).one()
            plate7443.plate_pointings[0].priority = 1

    def tearDown(self):
        """Sets default priorities for plates 8550 and 7443."""

        with session.begin(subtransactions=True):

            # Restores priorities
            plate8550 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 8550).one()
            plate8550.plate_pointings[0].priority = 5

            plate7443 = session.query(db.plateDB.Plate).filter(
                db.plateDB.Plate.plate_id == 7443).one()
            plate7443.plate_pointings[0].priority = 1

            # Restores the list of active pluggings
            activePluggings = session.query(db.plateDB.ActivePlugging).all()
            for actPlug in activePluggings:
                session.delete(actPlug)

            for originalActivePlugging in self.originalActivePluggings:
                pk = originalActivePlugging.pk
                plugging_pk = originalActivePlugging.plugging_pk
                session.add(db.plateDB.ActivePlugging(pk=pk,
                                                      plugging_pk=plugging_pk))

        # Restores offline carts
        config['offlineCarts'] = self.originalOfflineCarts

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

    def testCartOrder(self):
        """Tests the cart order returned by the Plugger."""

        # Sets a custom list of active pluggings.
        with session.begin(subtransactions=True):

            for actPlug in self.originalActivePluggings:
                session.delete(actPlug)

            session.add(db.plateDB.ActivePlugging(plugging_pk=68904, pk=2))
            session.add(db.plateDB.ActivePlugging(plugging_pk=70003, pk=1))
            session.add(db.plateDB.ActivePlugging(plugging_pk=70124, pk=3))

        # Calls the Plugger without dates
        plugger = Plugger(startDate=None, endDate=None)

        self.assertEqual(plugger.getASOutput()['cart_order'],
                         [9, 8, 7, 4, 5, 6, 3, 2, 1])

        # Restores list of original active plugging. For extra safety, this is
        # also done in tearDown.
        with session.begin(subtransactions=True):

            activePluggings = session.query(db.plateDB.ActivePlugging).all()
            for actPlug in activePluggings:
                session.delete(actPlug)

            for originalActivePlugging in self.originalActivePluggings:
                pk = originalActivePlugging.pk
                plugging_pk = originalActivePlugging.plugging_pk
                session.add(db.plateDB.ActivePlugging(pk=pk,
                                                      plugging_pk=plugging_pk))

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
                                                      2, 3, 1])]))

if __name__ == '__main__':
    unittest.main()
