#!/usr/bin/env python
# encoding: utf-8
"""
testSetArrangement.py

Created by José Sánchez-Gallego on 14 Jul 2014.
Licensed under a 3-clause BSD license.

Revision history:
    3 May 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from Totoro import TotoroDBConnection
from Totoro.dbclasses import fromPlateID
from Totoro.dbclasses.plate_utils import removeOrphanedSets
import unittest


db = TotoroDBConnection()
session = db.session


class TestSetArrangement(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Restores plate 7495."""

        setPK = [1, 1, 1, 2, 2, 2, 3, 4, 3, 4, 3, 4]

        with db.session.begin():

            for ii, expPK in enumerate(range(17, 29)):
                exp = db.session.query(db.mangaDB.Exposure).get(expPK)
                ss = db.session.query(db.mangaDB.Set).get(setPK[ii])
                if ss is None:
                    db.session.add(db.mangaDB.Set(pk=setPK[ii]))
                    db.session.flush()
                exp.set_pk = setPK[ii]
                exp.exposure_status_pk = 4
                db.session.flush()

            for sPK in setPK:
                ss = db.session.query(db.mangaDB.Set).get(sPK)
                ss.set_status_pk = 0

        removeOrphanedSets()

    @classmethod
    def tearDownClass(cls):
        """Similar to set up."""

        cls.setUpClass()

    def testExposureAssignment(self):
        """Tests if an exposure is assigned to a correct incomplete set."""

        with session.begin():
            exposure = session.query(db.mangaDB.Exposure).get(1348)
            exposure.set_pk = None

        plate = fromPlateID(8484, rearrangeIncomplete=False, force=True)

        assert len(plate.sets) == 6

        setExposurePK = [exp._mangaExposure.pk
                         for exp in plate.sets[4].totoroExposures]
        self.assertIn(1348, setExposurePK)

    def testIncompleteSetsRearrangement(self):
        """Tests whether the rearrang. of incomplete sets works properly."""

        # Removes all set assignment for plate 8551
        plate8551 = fromPlateID(8551)
        with session.begin():
            exposures = plate8551.getScienceExposures()
            for exp in exposures:
                setPK = exp.mangadbExposure[0].set_pk
                exp.mangadbExposure[0].set_pk = None
                if setPK is not None:
                    ss = session.query(db.mangaDB.Set).get(setPK)
                    session.delete(ss)

        # Reloads plate 8551
        plate8551 = fromPlateID(8551, force=True)

        # Checks new arrangement
        self.assertEqual(len(plate8551.sets), 5)

        # This is the expected assignment of exposures for each set
        correctSetExposures = [[198628, 198629],
                               [198624],
                               [198625, 198626, 198627],
                               [198623, 198621, 198622],
                               [198618, 198619, 198620]]
        for ii, ss in enumerate(plate8551.sets):
            setExposures = [exp.exposure_no
                            for exp in plate8551.sets[ii].totoroExposures]
            self.assertItemsEqual(setExposures, correctSetExposures[ii])

    def testRearrangementWithOverriddenSet(self):
        """Tests a rearrangement in a plate with an overriden set."""

        with db.session.begin():
            exp18 = db.session.query(db.mangaDB.Exposure).get(18)
            exp20 = db.session.query(db.mangaDB.Exposure).get(20)
            exp18.set_pk = 2
            exp20.set_pk = 1

            set1 = db.session.query(db.mangaDB.Set).get(1)
            set1.set_status_pk = 3

        plate = fromPlateID(7495)
        plate.rearrangeSets()

        correctSetExposures = [[17, 20, 19],
                               [18],
                               [23, 21, 22],
                               [24, 25],
                               [27, 28, 26]]

        for ii, ss in enumerate(plate.sets):
            setExposures = [exp._mangaExposure.pk
                            for exp in plate.sets[ii].totoroExposures]
            self.assertItemsEqual(setExposures, correctSetExposures[ii])


if __name__ == '__main__':
    unittest.main()
