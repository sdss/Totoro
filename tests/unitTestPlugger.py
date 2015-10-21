#!/usr/bin/env python
# encoding: utf-8
"""
unitTestPlugger.py

Created by José Sánchez-Gallego on 12 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    12 Aug 2014 J. Sánchez-Gallego
      Initial version
    20 Oct 2014
      Converted to testcase

"""

from __future__ import division
from __future__ import print_function
import unittest
from sdss.internal.manga.Totoro.scheduler import Plugger


class pluggerTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.plugger = Plugger(2456951.82708, 2456952.00208)

    def testSchedule(self):

        observingBlocks = self.plugger.observingBlocks

        self.assertEqual(len(observingBlocks), 1)
        self.assertAlmostEqual(observingBlocks[0]['JD0'], 2456951.82708)
        self.assertAlmostEqual(observingBlocks[0]['JD1'], 2456952.00208)

    def testPlates(self):

        platesAtAPO = self.plugger.plates
        plateIDs = [plate.plate_id for plate in platesAtAPO]

        self.assertGreater(len(platesAtAPO), 1)
        self.assertIn(7341, plateIDs)
        self.assertIn(8082, plateIDs)

    def testPlugger(self):

        pluggerSchedule = self.plugger.getOutput()
        print(pluggerSchedule.timeline.lstRange)


if __name__ == '__main__':
    unittest.main()
