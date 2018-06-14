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

from __future__ import division, print_function

import unittest
import warnings

import numpy as np

from Totoro.bin.overrideSet import main
from Totoro.db import getConnectionFull
from Totoro.dbclasses import Exposure, Set, fromPlateID
from Totoro.dbclasses.plate_utils import removeOrphanedSets


db, Session, plateDB, mangaDB = getConnectionFull('test')
session = db.Session()


class TestOverrideSet(unittest.TestCase):

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
                    session.flush()
                exp.set_pk = setPK[ii]
                exp.exposure_status_pk = 4
                session.flush()

            for sPK in setPK:
                ss = session.query(db.mangaDB.Set).get(sPK)
                ss.set_status_pk = 0

        plate8549_setPKs = [268, 269, 276, 280, 281]
        plate8549_setStatus = [1, None, 1, None, 1]
        with session.begin():
            for ii, setPK in enumerate(plate8549_setPKs):
                ss = session.query(db.mangaDB.Set).get(setPK)
                if ss is None:
                    session.add(db.mangaDB.Set(pk=setPK, set_status_pk=plate8549_setStatus[ii]))
                    session.flush()
                    continue
                else:
                    ss.set_status_pk = plate8549_setStatus[ii]
                    session.flush()

                for exp in ss.exposures:
                    exp.set_pk = None

        plate8549_exposures = [
            84465, 84588, 84667, 84468, 84462, 84664, 84670, 84658, 84661, 84591, 84672, 84675
        ]
        plate8549_exposure_set = [268, 268, 268, 269, 276, 276, 276, 280, 280, 281, 281, 281]
        with session.begin():
            for ii, expPK in enumerate(plate8549_exposures):
                exp = session.query(db.plateDB.Exposure).get(expPK)
                if expPK == 84468:
                    exp.mangadbExposure[0].exposure_status_pk = None
                else:
                    exp.mangadbExposure[0].exposure_status_pk = 4
                exp.mangadbExposure[0].set_pk = plate8549_exposure_set[ii]

            ss = session.query(db.mangaDB.Set).get(312)
            if ss is not None:
                session.delete(ss)

        removeOrphanedSets()

    def tearDown(self):
        """Similar to set up."""

        self.setUpClass()

    def testOverrideGoodOneSet(self):
        """Test overriding all the exposures in a set to Override Good."""

        exps = [177773, 177774, 177775]
        main(argv=['-v', 'good'] + list(map(str, exps)))

        for exp in exps:
            totExp = Exposure(exp, format='exposure_no', parent='plateDB')
            self.assertEqual(totExp._mangaExposure.set_pk, 1)
            self.assertEqual(totExp._mangaExposure.status.label, 'Override Good')
            self.assertEqual(totExp._mangaExposure.set.status.label, 'Override Good')

        # Checks that the set is overriden good
        with session.begin():
            ss = session.query(db.mangaDB.Set).get(1)

        self.assertIsNotNone(ss)
        self.assertEqual(ss.status.label, 'Override Good')

    def testOverrideGoodSeveralSet(self):
        """Test overriding exposures from different sets into a good set."""

        exps = [177773, 177774, 177778]
        main(argv=['-v', 'good'] + list(map(str, exps)))

        setPKs = []
        for exp in exps:
            totExp = Exposure(exp, format='exposure_no', parent='plateDB')
            setPK = totExp._mangaExposure.set_pk
            self.assertIsNotNone(setPK)
            setPKs.append(setPK)
            self.assertEqual(totExp._mangaExposure.status.label, 'Override Good')
            self.assertEqual(totExp._mangaExposure.set.status.label, 'Override Good')

        self.assertEqual(len(np.unique(setPKs)), 1)

        # Checks that the set is overriden good
        with session.begin():
            ss = session.query(db.mangaDB.Set).get(setPKs[0])

        self.assertIsNotNone(ss)
        self.assertEqual(ss.status.label, 'Override Good')

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
        main(argv=['-v', 'bad'] + list(map(str, exps)))

        for exp in exps:
            totExp = Exposure(exp, format='exposure_no', parent='plateDB')
            self.assertEqual(totExp._mangaExposure.set_pk, 1)
            self.assertEqual(totExp._mangaExposure.status.label, 'Override Bad')
            self.assertEqual(totExp._mangaExposure.set.status.label, 'Override Bad')

        # Checks that the set is overriden good
        with session.begin():
            ss = session.query(db.mangaDB.Set).get(1)

        self.assertIsNotNone(ss)
        self.assertEqual(ss.status.label, 'Override Bad')

        plate = fromPlateID(7495)
        self.assertEqual(len(plate.sets), 4)

        # Checks that the plate completion takes into account that there is one
        # fewer valid set.
        self.assertLessEqual(plate.getPlateCompletion(), 1.5)

    def testOverrideBadSeveralSet(self):
        """Test overriding exposures from different sets into a bad set."""

        exps = [177773, 177774, 177778]

        with warnings.catch_warnings(record=True) as ww:

            main(argv=['-v', 'bad'] + list(map(str, exps)))

            warnMessages = '\n'.join([str(ww[ii].message) for ii in range(len(ww))])

            self.assertIn('plate completion has changed from 1.85 to 0.91.', warnMessages)

        setPKs = []
        for exp in exps:
            totExp = Exposure(exp, format='exposure_no', parent='plateDB')
            setPK = totExp._mangaExposure.set_pk
            self.assertIsNotNone(setPK)
            setPKs.append(setPK)
            self.assertEqual(totExp._mangaExposure.status.label, 'Override Bad')
            self.assertEqual(totExp._mangaExposure.set.status.label, 'Override Bad')

        self.assertEqual(len(np.unique(setPKs)), 1)

        # Checks that the set is overriden bad
        with session.begin():
            ss = session.query(db.mangaDB.Set).get(setPKs[0])

        self.assertIsNotNone(ss)
        self.assertEqual(ss.status.label, 'Override Bad')

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
        overridenSetPK = main(argv=['-v', 'bad'] + list(map(str, exps)))

        # Checks that the new overridden set is 297
        with session.begin():
            ss = session.query(db.mangaDB.Set).get(overridenSetPK)

        self.assertIsNotNone(ss)
        self.assertEqual(ss.status.label, 'Override Bad')

        # Now we remove the set
        main(argv=['-v', 'remove', str(overridenSetPK)])

        # Check that set_pk 297 has been removed
        with session.begin():
            ss = session.query(db.mangaDB.Set).get(overridenSetPK)

        self.assertIsNone(ss)

        # Checks that all exposures don't have exposure status or set_pk
        with session.begin():
            for exp in exps:
                ee = session.query(
                    db.plateDB.Exposure).filter(db.plateDB.Exposure.exposure_no == exp).one()
                mangaDBexp = ee.mangadbExposure[0]
                self.assertIsNone(mangaDBexp.set_pk)
                self.assertIsNone(mangaDBexp.exposure_status_pk)

        # Now we repeat the test but using --reload. The plate should be left
        # in the original state but with sets 1 and 2 exchanged.
        overridenSetPK = main(argv=['-v', 'bad'] + list(map(str, exps)))

        # Now we remove the set
        main(argv=['-v', 'remove', '--reload', str(overridenSetPK)])

        # Checks that all exposures don't have exposure status or set_pk
        expSetPKs = [1, 1, 2]
        with session.begin():
            for ii, exp in enumerate(exps):
                ee = session.query(
                    db.plateDB.Exposure).filter(db.plateDB.Exposure.exposure_no == exp).one()
                mangaDBexp = ee.mangadbExposure[0]
                self.assertEqual(mangaDBexp.set_pk, expSetPKs[ii])
                self.assertEqual(mangaDBexp.exposure_status_pk, 4)

        setStatus = [1, None]
        with session.begin():
            for ii, setPK in enumerate([1, 2]):
                ss = session.query(db.mangaDB.Set).get(setPK)
                self.assertEqual(ss.set_status_pk, setStatus[ii])

    def testInfo(self):
        """Test getting information from a set."""

        # Checks a fake set
        exps = [177773, 177774, 177778]
        (status, code, statusMock, codeMock) = main(argv=['-v', 'info'] + list(map(str, exps)))

        self.assertEqual(status, 'Bad')
        self.assertEqual(code, 2)
        self.assertIsNone(statusMock)
        self.assertIsNone(codeMock)

        # Now we check a real good one
        exps = [177773, 177774, 177775]
        (status, code, statusMock, codeMock) = main(argv=['-v', 'info'] + list(map(str, exps)))

        self.assertEqual(status, 'Excellent')
        self.assertEqual(code, 10)
        self.assertEqual(statusMock, 'Good')
        self.assertEqual(codeMock, 0)

        # Lets override the first example as bad
        exps = [177773, 177774, 177778]
        main(argv=['-v', 'bad'] + list(map(str, exps)))
        (status, code, statusMock, codeMock) = main(argv=['-v', 'info'] + list(map(str, exps)))

        self.assertEqual(status, 'Override Bad')
        self.assertEqual(code, 10)
        self.assertEqual(statusMock, 'Bad')
        self.assertEqual(codeMock, 2)

        # For the previous test the status of the exposures has been internally
        # changed to None. That change should not be recorded in the DB. Let's
        # check.
        for exp in exps:
            totExp = Exposure(exp, format='exposure_no')
            self.assertEqual(totExp._mangaExposure.status.label, 'Override Bad')

    def testEmptySet(self):
        """Tests when one of the original sets in empty after overriding."""

        exps = [198371, 198447, 198258]
        main(argv=['-v', 'good'] + list(map(str, exps)))

        with session.begin():
            exp198371 = session.query(
                db.plateDB.Exposure).filter(db.plateDB.Exposure.exposure_no == 198371).one()
            ss = exp198371.mangadbExposure[0].set

        self.assertIsNotNone(ss)
        self.assertEqual(ss.status.label, 'Override Good')
        self.assertItemsEqual([exp.platedbExposure.exposure_no for exp in ss.exposures], exps)


if __name__ == '__main__':
    unittest.main()
