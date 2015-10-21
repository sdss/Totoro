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
from sdss.internal.manga.Totoro.dbclasses import Plate, Exposure, Set
import numpy as np


class testSets():

    def testSetLoad(self):
        """Tests set loading."""

        plateID = 7815
        plate = Plate(plateID, format='plate_id')

        assert(len(plate.sets) > 0)

    def testSetCoordinates(self):
        """Tests set coordinates."""

        plateID = 7815
        plate = Plate(plateID, format='plate_id')

        np.testing.assert_almost_equal(plate.ra, plate.sets[0].ra)
        np.testing.assert_almost_equal(plate.dec, plate.sets[0].dec)

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

        assert nSet is not None, 'could not find set to perform test'

        assert plate8484.sets[nSet].getStatus()[0] == 'Unplugged'

        # Tests creating a new mock exposure and adding it to an unplugged set.
        # The status of the new set should be bad because the pluggings are
        # different.
        newExp = Exposure.createMockExposure(
            startTime=2457135.9304166664, expTime=900,
            ra=plate8484.ra, dec=plate8484.dec, ditherPosition='S')
        exposures = plate8484.sets[nSet].totoroExposures + [newExp]
        newSet = Set.fromExposures(exposures)

        assert newSet.getStatus()[0] == 'Bad'

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

        assert newSet2.getStatus()[0] == 'Good'
