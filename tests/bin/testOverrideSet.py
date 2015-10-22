#!/usr/bin/env python
# encoding: utf-8
"""
testOverrideSet.py

Created by José Sánchez-Gallego on 27 Sep 2015.
Licensed under a 3-clause BSD license.

Revision history:
    27 Sep 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import unittest
import numpy as np
import warnings

from sdss.internal.manga.Totoro import TotoroDBConnection
from sdss.internal.manga.Totoro.dbclasses.plate_utils import removeOrphanedSets
from sdss.internal.manga.Totoro.bin.overrideSet import main
from sdss.internal.manga.Totoro.dbclasses import Exposure, Set, fromPlateID


db = TotoroDBConnection()


class TestOverrideSet(unittest.TestCase):

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

    def tearDown(self):
        """Similar to set up."""

        self.setUpClass()

    def testOverrideGoodOneSet(self):
        """Test overriding all the exposures in a set to Override Good."""

        exps = [177773, 177774, 177775]
        main(argv=['-v', 'good'] + map(str, exps))

        for exp in exps:
            totExp = Exposure(exp, format='exposure_no', parent='plateDB')
            self.assertEqual(totExp._mangaExposure.set_pk, 1)
            self.assertEqual(totExp._mangaExposure.status.label,
                             'Override Good')
            self.assertEqual(totExp._mangaExposure.set.status.label,
                             'Override Good')

    def testOverrideGoodSeveralSet(self):
        """Test overriding exposures from different sets into a good set."""

        exps = [177773, 177774, 177778]
        main(argv=['-v', 'good'] + map(str, exps))

        setPKs = []
        for exp in exps:
            totExp = Exposure(exp, format='exposure_no', parent='plateDB')
            setPK = totExp._mangaExposure.set_pk
            self.assertIsNotNone(setPK)
            setPKs.append(setPK)
            self.assertEqual(totExp._mangaExposure.status.label,
                             'Override Good')
            self.assertEqual(totExp._mangaExposure.set.status.label,
                             'Override Good')

        self.assertEqual(len(np.unique(setPKs)), 1)

        for setPK in [1, 2, 3, 4]:
            ss = Set(setPK, format='pk')
            status = ss.getStatus()[0]
            if setPK in [1, 2]:
                self.assertEqual(status, 'Unplugged')
            else:
                self.assertEqual(status, 'Good')
                self.assertEqual(ss.status.label, 'Excellent')

        plate = fromPlateID(7495)
        self.assertEqual(len(plate.sets), 5)

    def testOverrideBadOneSet(self):
        """Test overriding all the exposures in a set to Override Bad."""

        exps = [177773, 177774, 177775]
        main(argv=['-v', 'bad'] + map(str, exps))

        for exp in exps:
            totExp = Exposure(exp, format='exposure_no', parent='plateDB')
            self.assertEqual(totExp._mangaExposure.set_pk, 1)
            self.assertEqual(totExp._mangaExposure.status.label,
                             'Override Bad')
            self.assertEqual(totExp._mangaExposure.set.status.label,
                             'Override Bad')

        plate = fromPlateID(7495)
        self.assertEqual(len(plate.sets), 4)

        # Checks that the plate completion takes into account that there is one
        # fewer valid set.
        self.assertLessEqual(plate.getPlateCompletion(), 1.5)

    def testOverrideBadSeveralSet(self):
        """Test overriding exposures from different sets into a bad set."""

        exps = [177773, 177774, 177778]

        with warnings.catch_warnings(record=True) as ww:

            main(argv=['-v', 'bad'] + map(str, exps))

            warnMessages = '\n'.join([str(ww[ii].message)
                                      for ii in range(len(ww))])

            self.assertIn('plate completion has changed from 1.85 to 0.91.',
                          warnMessages)

        setPKs = []
        for exp in exps:
            totExp = Exposure(exp, format='exposure_no', parent='plateDB')
            setPK = totExp._mangaExposure.set_pk
            self.assertIsNotNone(setPK)
            setPKs.append(setPK)
            self.assertEqual(totExp._mangaExposure.status.label,
                             'Override Bad')
            self.assertEqual(totExp._mangaExposure.set.status.label,
                             'Override Bad')

        self.assertEqual(len(np.unique(setPKs)), 1)

        for setPK in [1, 2, 3, 4]:
            ss = Set(setPK, format='pk')
            status = ss.getStatus()[0]
            if setPK in [1, 2]:
                self.assertEqual(status, 'Unplugged')
            else:
                self.assertEqual(status, 'Good')
                self.assertEqual(ss.status.label, 'Excellent')

        plate = fromPlateID(7495)
        self.assertEqual(len(plate.sets), 5)

        self.assertLessEqual(plate.getPlateCompletion(), 1.)

    def testRemoveSet(self):
        """Test removing overridden set."""

        # We override a set as bad
        exps = [177773, 177774, 177778]
        main(argv=['-v', 'bad'] + map(str, exps))

        # Checks that the new overridden set is 297
        with db.session.begin():
            ss = db.session.query(db.mangaDB.Set).get(297)

        self.assertIsNotNone(ss)
        self.assertEqual(ss.status.label, 'Override Bad')

        # Now we remove the set
        main(argv=['-v', 'remove', '297'])

        # Check that set_pk 297 has been removed
        with db.session.begin():
            ss = db.session.query(db.mangaDB.Set).get(297)

        self.assertIsNone(ss)

        # Checks that all exposures don't have exposure status or set_pk
        with db.session.begin():
            for exp in exps:
                ee = db.session.query(db.plateDB.Exposure).filter(
                    db.plateDB.Exposure.exposure_no == exp).one()
                mangaDBexp = ee.mangadbExposure[0]
                self.assertIsNone(mangaDBexp.set_pk)
                self.assertIsNone(mangaDBexp.exposure_status_pk)

        # Now we repeat the test but using --reload. The plate should be left
        # in the original state but with sets 1 and 2 exchanged.
        main(argv=['-v', 'bad'] + map(str, exps))

        # Now we remove the set
        main(argv=['-v', 'remove', '--reload', '297'])

        # Checks that all exposures don't have exposure status or set_pk
        expSetPKs = [2, 2, 1]
        with db.session.begin():
            for ii, exp in enumerate(exps):
                ee = db.session.query(db.plateDB.Exposure).filter(
                    db.plateDB.Exposure.exposure_no == exp).one()
                mangaDBexp = ee.mangadbExposure[0]
                self.assertEqual(mangaDBexp.set_pk, expSetPKs[ii])
                self.assertEqual(mangaDBexp.exposure_status_pk, 4)

        with db.session.begin():
            for setPK in [1, 2]:
                ss = db.session.query(db.mangaDB.Set).get(setPK)
                self.assertEqual(ss.status.label, 'Good')


if __name__ == '__main__':
    unittest.main()
