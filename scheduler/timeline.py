#!/usr/bin/env python
# encoding: utf-8
"""
timeline.py

Created by José Sánchez-Gallego on 1 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    1 Aug 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from astropy.time import Time
from sdss.internal.manga.Totoro import log, TotoroDBConnection
from sdss.internal.manga.Totoro import utils
from sdss.internal.manga.Totoro import logic
import numpy as np


db = TotoroDBConnection()
session = db.Session()


class Timelines(list):
    """A list of Timelines."""

    def __init__(self, observingBlocks, mode='planner', **kwargs):

        self.mode = mode
        initList = []
        for row in observingBlocks:
            initList.append(Timeline(row['JD0'], row['JD1']))

        list.__init__(self, initList)

    def schedule(self, plates, mode=None, **kwargs):
        """Schedules all timelines."""

        if mode is None:
            mode = self.mode

        for timeline in self:
            timeline.schedule(plates, mode=mode, **kwargs)


class Timeline(object):
    """A timeline object that can be used to schedule one observing night."""

    def __init__(self, startDate, endDate, **kwargs):

        self.startDate = startDate
        self.endDate = endDate
        self.plates = []

        self.unallocatedPlateWindow = np.array([self.startDate, self.endDate])
        self.unallocatedExps = self.unallocatedPlateWindow.copy()

    def allocateJDs(self, plates=None, **kwargs):
        """Allocates observed JDs for a list of plates."""

        plates = self.plates if plates is None else plates

        nominalMJD = Time((self.startDate+self.endDate)/2., format='jd',
                          scale='tai').mjd

        nExp = 0  # number of new exposures

        for plate in plates:
            for exp in plate.getTotoroExposures():
                if not exp.isValid:
                    continue
                jdRange = exp.getJD()
                if (jdRange[1] >= self.startDate and
                        jdRange[0] <= self.endDate):
                    self.unallocatedExps = utils.removeInterval(
                        self.unallocatedExps, jdRange, wrapAt=None)
                    nExp += 1

            dateRangePlate = plate.getUTRange(mjd=nominalMJD,
                                              returnType='date')
            jdRangePlate = [dateRangePlate[0].jd, dateRangePlate[1].jd]

            self.unallocatedPlateWindow = utils.removeInterval(
                self.unallocatedPlateWindow, jdRangePlate, wrapAt=None)

            return nExp

    @property
    def remainingTime(self):
        """Returns the amount of unallocated time, in hours."""
        unallocatedExps = np.atleast_2d(self.unallocatedExps)
        if unallocatedExps.size == 0:
            return 0
        unallocatedTime = 0.0
        for interval in unallocatedExps:
            unallocatedTime += (interval[1]-interval[0])
        return unallocatedTime * 24.

    def schedule(self, plates, mode='plugger', allowComplete=False,
                 showUnobservedTimes=True, **kwargs):
        """Schedules a list of plates in the LST ranges not yet observed in the
        timeline."""

        prioritisePlugged = True if mode == 'plugger' else False

        log.debug('scheduling LST range {0} using {1} plates, mode={2}'
                  .format(self.unallocatedExps.tolist(), len(plates), mode))

        if not allowComplete:
            plates = [plate for plate in plates if not plate.isComplete]

        while self.remainingTime > 0:

            optimalPlate = logic.getOptimalPlate(
                plates, self.unallocatedExps,
                prioritisePlugged=prioritisePlugged, mode=mode)

            if optimalPlate is None:
                break
            else:
                self.plates.append(optimalPlate)
                nExp = self.allocateJDs(plates=[optimalPlate])
                log.info('...... found optimal plate: plate_id={0}, '
                         'manga_tiledid={1} ({2} new exposures)'
                         .format(optimalPlate.plate_id,
                                 optimalPlate.getMangaTileID(),
                                 nExp))

        if showUnobservedTimes:
            if self.remainingTime == 0:
                return True
            else:
                log.info('... unobserved times: {0}'.format(
                         str(self.unallocatedExps)))
                return False
