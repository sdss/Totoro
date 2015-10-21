#!/usr/bin/env python
# encoding: utf-8
"""
unitTestSets.py

Created by José Sánchez-Gallego on 27 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    27 Aug 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import unittest
from Totoro.dbclasses import Plate


class setTestCase(unittest.TestCase):

    def testSetLoad(self):

        plateID = 7815
        plate = Plate(plateID, format='plate_id', sets=True)

        self.assertGreater(len(plate.sets), 0)

    def testSetCoordinates(self):

        plateID = 7815
        plate = Plate(plateID, format='plate_id', sets=True)

        self.assertAlmostEqual(plate.ra, plate.sets[0].ra)
        self.assertAlmostEqual(plate.dec, plate.sets[0].dec)


if __name__ == '__main__':
    unittest.main()
