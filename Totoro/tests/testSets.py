#!/usr/bin/env python
# encoding: utf-8
"""
testSets.py

Created by José Sánchez-Gallego on 27 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    27 Aug 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from builtins import range
from Totoro.dbclasses import Plate, Exposure, Set
from Totoro.db import getConnection
from Totoro.dbclasses.plate_utils import removeOrphanedSets
import numpy as np
import unittest

db = getConnection('test')
session = db.Session()


class TestSets(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Restores plate 7495."""

        setPK = [1, 1, 1, 2, 2, 2, 3, 4, 3, 4, 3, 4]

        with session.begin():

            for ii, expPK in enumerate(range(17, 29)):
                exp = session.query(db.mangaDB.Exposure).get(expPK)
                ss = session.query(db.mangaDB.Set).get(setPK[ii])
                if ss is None:
                    session.add(db.mangaDB.Set(pk=setPK[ii]))
                exp.set_pk = setPK[ii]
                exp.exposure_status_pk = 4

            sn2 = session.query(db.mangaDB.SN2Values).get(17)
            sn2.b1_sn2 = 3.28888
            sn2.b2_sn2 = 2.90681
            sn2.r1_sn2 = 6.01449
            sn2.r2_sn2 = 6.55439

        removeOrphanedSets()

    @classmethod
    def tearDownClass(cls):
        """Similar to set up."""

        cls.setUpClass()

    def testSetLoad(self):
        """Tests set loading."""

        plateID = 7815
        plate = Plate(plateID, format='plate_id')

        self.assertGreater(len(plate.sets), 0)

    def testSetCoordinates(self):
        """Tests set coordinates."""

        plateID = 7815
        plate = Plate(plateID, format='plate_id')

        self.assertAlmostEqual(plate.ra, plate.sets[0].ra)
        self.assertAlmostEqual(plate.dec, plate.sets[0].dec)

    def testUnplugged(self):
        """Tests behaviour of unplugged sets and mock sets."""

        plate8484 = Plate(8484, format='plate_id')

        # Finds the unplugged set with two exposures. This set may change
        # depending on other tests run before
        nSet = None
        for nn, ss in enumerate(plate8484.sets):
            if (ss.getStatus()[0] == 'Unplugged' and
                    len(ss.totoroExposures) == 2):
                nSet = nn
                break

        self.assertNotEqual(nSet, None)

        self.assertEqual(plate8484.sets[nSet].getStatus()[0], 'Unplugged')

        # Tests creating a new mock exposure and adding it to an unplugged set.
        # The status of the new set should be bad because the pluggings are
        # different.
        newExp = Exposure.createMockExposure(
            startTime=2457135.9304166664, expTime=900,
            ra=plate8484.ra, dec=plate8484.dec, ditherPosition='S')
        exposures = plate8484.sets[nSet].totoroExposures + [newExp]
        newSet = Set.fromExposures(exposures)

        self.assertEqual(newSet.getStatus()[0], 'Bad')

        # Repeats the same test but now the mock exposure is created with the
        # plugging of one of the exposures in the unplugged set.
        plugging = plate8484.sets[nSet].totoroExposures[0].getPlugging()

        newExp2 = Exposure.createMockExposure(
            startTime=2457135.9304166664, expTime=900,
            ra=plate8484.ra, dec=plate8484.dec, ditherPosition='S',
            plugging=plugging)

        # Tricks SN2 values to avoid the set failing because of SN2 uniformity
        newExp2._sn2Array = np.array([1.5, 1.5, 3, 3])

        # Creates new mock set
        exposures = plate8484.sets[nSet].totoroExposures + [newExp2]
        newSet2 = Set.fromExposures(exposures)

        self.assertEqual(newSet2.getStatus()[0], 'Good')

        # Now let's use a plate with active pluggings to check the status
        # of new mock exposures.

        plate8486 = Plate(8486, format='plate_id')

        self.assertEqual(len(plate8486.sets), 2)
        for ss in plate8486.sets:
            self.assertEqual(ss.getStatus()[0], 'Good')

        plate8486.addMockExposure(startTime=2457137.9535648148)
        self.assertEqual(len(plate8486.sets), 3)
        statuses = [ss.getStatus()[0] for ss in plate8486.sets]
        self.assertEqual(statuses.count('Incomplete'), 1)

        plate8486.addMockExposure(startTime=2457137.9535648148,
                                  plugging=plate8486.getActivePlugging())
        self.assertEqual(len(plate8486.sets), 4)
        statuses = [ss.getStatus()[0] for ss in plate8486.sets]
        self.assertEqual(statuses.count('Incomplete'), 2)

        plate8486.addMockExposure(startTime=2457137.9535648148,
                                  plugging=plate8486.getActivePlugging())
        self.assertEqual(len(plate8486.sets), 4)
        statuses = [ss.getStatus()[0] for ss in plate8486.sets]
        self.assertEqual(statuses.count('Incomplete'), 2)

    def testIncompleteExposure(self):
        """Checks if an incompletely reduced exposure behaves properly."""

        # Checks exposure
        exp = Exposure(17, parent='mangadb')
        self.assertIsNot(exp._mangaExposure.set_pk, None)

        # Modifies the exposure
        with session.begin():
            exp = session.query(db.mangaDB.Exposure).get(17)
            exp.set_pk = None
            exp.sn2values[0].b1_sn2 = None
            exp.sn2values[0].b2_sn2 = None
            exp.sn2values[0].r1_sn2 = None
            exp.sn2values[0].r2_sn2 = None
            exp.exposure_status_pk = None

        exp = Exposure(17, parent='mangadb')
        self.assertIsNone(exp._mangaExposure.set_pk)

        # Loads the plate
        plate = Plate(7495, format='plate_id', force=True)
        for exp in plate.getScienceExposures():
            if exp.mangadbExposure[0].pk == 17:
                self.assertIsNone(exp.mangadbExposure[0].set_pk)
                break

        # Changes exposure to incomplete but valid
        with session.begin():
            exp = session.query(db.mangaDB.Exposure).get(17)
            exp.sn2values[0].b1_sn2 = 3.28888
            exp.sn2values[0].b2_sn2 = 2.90681
            exp.sn2values[0].r1_sn2 = 6.01449
            exp.sn2values[0].r2_sn2 = None

        plate = Plate(7495, format='plate_id', force=True)
        for exp in plate.getTotoroExposures():
            if exp._mangaExposure.pk == 17:
                # self.assertEqual(exp.mangadbExposure[0].set_pk, 1)
                self.assertIsNone(exp.mangadbExposure[0].status)
                break

        # Now we break the exposure again
        with session.begin():
            exp = session.query(db.mangaDB.Exposure).get(17)
            exp.sn2values[0].b1_sn2 = 0.0
            exp.sn2values[0].b2_sn2 = 0.0
            exp.sn2values[0].r1_sn2 = 6.01449
            exp.sn2values[0].r2_sn2 = 6.55439

        # plate = Plate(7495, format='plate_id', force=True)
        # for exp in plate.getTotoroExposures():
        #     if exp._mangaExposure.pk == 17:
        #         self.assertIsNone(exp.mangadbExposure[0].set_pk)
        #         self.assertEqual(exp.mangadbExposure[0].status.pk, 5)
        #         break

if __name__ == '__main__':
    unittest.main()
