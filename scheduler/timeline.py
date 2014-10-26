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
import warnings
from astropy.time import Time
from sdss.internal.manga.Totoro import config, log, site, TotoroDBConnection
from sdss.internal.manga.Totoro.exceptions import TotoroUserWarning
from sdss.internal.manga.Totoro import utils
from sdss.internal.manga.Totoro import logic
import numpy as np


db = TotoroDBConnection()
session = db.Session()


class Timelines(list):
    """A list of Timelines."""

    def __init__(self, observingBlocks, plates=None, mode='planner', **kwargs):

        self.mode = mode
        initList = []
        for row in observingBlocks:
            initList.append(Timeline(row['JD0'], row['JD1'], plates=plates))

        list.__init__(self, initList)

    def schedule(self, mode=None, **kwargs):
        """Schedules all timelines."""

        if mode is None:
            mode = self.mode

        plates = self[0]._plates

        for timeline in self:
            timeline._plates = plates
            timeline.schedule(mode=mode, **kwargs)
            plates = timeline._plates


class Timeline(object):
    """A timeline object that can be used to schedule one observing night."""

    def __init__(self, startTime, endTime, **kwargs):

        self.startTime = startTime
        self.endTime = endTime
        self._plates = []

        self.unallocatedPlateWindow = np.array([self.startTime, self.endTime])
        self.unallocatedExps = self.unallocatedPlateWindow.copy()

    def allocateJDs(self, plates=None, **kwargs):
        """Allocates observed JDs for a list of plates."""

        plates = self._plates if plates is None else plates

        nominalMJD = int(Time((self.startTime+self.endTime)/2., format='jd',
                              scale='tai').mjd)

        for plate in plates:
            for exp in plate.getTotoroExposures():
                if not exp.isValid:
                    continue
                jdRange = exp.getJD()
                if (jdRange[1] >= self.startTime and
                        jdRange[0] <= self.endTime):
                    self.unallocatedExps = utils.removeInterval(
                        self.unallocatedExps, jdRange, wrapAt=None)

            dateRangePlate = plate.getUTRange(mjd=nominalMJD,
                                              returnType='date')
            jdRangePlate = [dateRangePlate[0].jd, dateRangePlate[1].jd]
            self.unallocatedPlateWindow = utils.removeInterval(
                self.unallocatedPlateWindow, jdRangePlate, wrapAt=None)

    def schedule(self, plates, mode='plugger', force=False,
                 **kwargs):
        """Schedules a list of plates in the LST ranges not yet observed in the
        timeline."""

        prioritisePlugged = True if mode == 'plugger' else False

        log.debug('scheduling LST range {0} using {1} plates, '
                  'mode={2}, force={3}'
                  .format(self.unallocatedExps.tolist(),
                          len(plates), mode, force))

        if not force:
            plates = [plate for plate in plates if not plate.isComplete]

        while self.unallocatedExps.size > 0 and len(plates) > 0:

            optimalPlate = logic.getOptimalPlate(
                plates, self.unallocatedExps,
                prioritisePlugged=prioritisePlugged)

            if optimalPlate is None:
                break
            else:
                self._plates.append(optimalPlate)
                self.allocateJDs(plates=[optimalPlate])
                plates.remove(optimalPlate)

        if len(plates) > 0 and force:
            self.allocateJDs(plates)
            self._plates += plates

        if self.unallocatedExps.size == 0:
            return True
        else:
            return False

    # def schedule(self, mode='planner', **kwargs):

    #     log.info('scheduling timeline with JD0={0:.4f}, JD1={1:.4f}'
    #              .format(self.startTime, self.endTime))

    #     currentTime = self.startTime
    #     remainingTime = utils.JDdiff(currentTime, self.endTime)

    #     # expTime = (config['exposure']['exposureTime'] /
    #     #            config[mode]['efficiency'])

    #     nCart = 1
    #     totalCarts = config[mode]['nCarts']

    #     while remainingTime >= 0. and nCart <= totalCarts:

    #         log.info('... simulating plates for {0:.5f}'.format(currentTime))

    #         optimalPlate = logic.getOptimalPlate(
    #             self._plates, currentTime, self.endTime,
    #             mode=mode, **kwargs)

    #         if optimalPlate is None:

    #             warnings.warn('no valid plates found at JD={0:.4f} '
    #                           '(timeline ends at JD={1:.4f})'.format(
    #                               currentTime, self.endTime),
    #                           TotoroUserWarning)
    #             currentTime += config['exposure']['exposureTime'] / 86400.

    #         else:

    #             nCart += 1
    #             currentTime = optimalPlate.getLastExposure().getJD()[1]

    #         remainingTime = utils.JDdiff(currentTime, self.endTime)

    #     if nCart > totalCarts:
    #         warnings.warn('run out of cart with {0} seconds remaining'
    #                       .format(remainingTime), TotoroUserWarning)

    #     return
