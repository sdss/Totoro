#!/usr/bin/env python
# encoding: utf-8
"""
testPlates.py

Created by José Sánchez-Gallego on 27 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    27 Aug 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import unittest
from sdss.internal.manga.Totoro.dbclasses import Plate, fromPlateID
from sdss.internal.manga.Totoro import TotoroDBConnection
from sdss.internal.manga.Totoro import exceptions


class plateTestCase(unittest.TestCase):

    def testPlateLoadWithPK(self):
        """Tests loading a plate from pk."""

        platePK = 11147
        plate = Plate(platePK)

        self.assertEqual(plate.plate_id, 7990)

    def testPlateLoadWithPlateID(self):
        """Tests loading a plate from plateID."""

        plateID = 7815
        plate1 = Plate(plateID, format='plate_id')
        plate2 = fromPlateID(plateID)

        self.assertEqual(plate1.comment, '060-25_MGA')
        self.assertEqual(plate2.location_id, 3760)

    def testPlateCoordinates(self):
        """Tests plate coordinates."""

        plateID = 7815
        plate = fromPlateID(plateID, sets=False)

        self.assertAlmostEqual(plate.coords[0], 317.95449707, places=4)
        self.assertAlmostEqual(plate.coords[1], 10.1960287094, places=4)

    def testIsPlugged(self):
        """Tests if plate is plugged."""

        plate = fromPlateID(8486)
        self.assertEqual(plate.isPlugged, True)

    def testMangadbExposures(self):
        """Tests accessing mangaDB exposures."""

        plate = fromPlateID(7815, sets=False)
        self.assertEqual(len(plate.getMangadbExposures()), 18)
        self.assertEqual(len(plate.getMangadbExposures()),
                         len(plate.getScienceExposures()))

    def testSubtransactions(self):
        """Fails if trying to load a plate from within a subtransaction."""

        db = TotoroDBConnection()
        with self.assertRaises(exceptions.TotoroError):
            with db.session.begin():
                fromPlateID(7815)


if __name__ == '__main__':
    unittest.main()
