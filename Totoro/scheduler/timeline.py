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
from sdss.internal.manga.Totoro import log, config
from sdss.internal.manga.Totoro import utils
from sdss.internal.manga.Totoro.scheduler import scheduler_utils as logic
from sdss.internal.manga.Totoro.core.colourPrint import _color_text
from sdss.internal.manga.Totoro.exceptions import TotoroUserWarning
import numpy as np
import warnings

expTimeJD = config['exposure']['exposureTime'] / 86400.


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
        self.scheduled = []

        self.unallocatedRange = np.array([self.startDate, self.endDate])

    def allocateJDs(self, exposures, **kwargs):
        """Allocates observed JDs for a list of plates."""

        for exp in exposures:
            expJD = exp.getJD()
            self.unallocatedRange = utils.removeInterval(
                self.unallocatedRange, expJD, wrapAt=None)

        return

    def calculatePlateCompletion(self, plate, rejectExposures=[],
                                 useMock=True):
        """Calculates plate completion after rejecting exposures."""

        from sdss.internal.manga.Totoro.dbclasses import Plate, Set

        if len(rejectExposures) == 0:
            return plate.getPlateCompletion(useMock=useMock)

        mockPlate = Plate.createMockPlate(ra=plate.ra, dec=plate.dec)
        for ss in plate.sets:
            newSetExps = []
            for exp in ss.totoroExposures:
                if exp not in rejectExposures:
                    newSetExps.append(exp)
            if len(newSetExps) > 0:
                mockPlate.sets.append(Set.fromExposures(newSetExps))

        return mockPlate.getPlateCompletion(useMock=useMock)

    @property
    def remainingTime(self):
        """Returns the amount of unallocated time, in hours."""
        unallocatedRange = np.atleast_2d(self.unallocatedRange)
        if unallocatedRange.size == 0:
            return 0
        unallocatedTime = 0.0
        for interval in unallocatedRange:
            unallocatedTime += (interval[1]-interval[0])
        return unallocatedTime * 24.

    def schedule(self, plates, mode='plugger', showUnobservedTimes=True,
                 useDateAtAPO=True, **kwargs):
        """Schedules a list of plates in the LST ranges not yet observed in the
        timeline."""

        from sdss.internal.manga.Totoro.dbclasses import Field

        mode = mode.lower()

        log.debug('scheduling LST range {0} using {1} plates, mode={2}'
                  .format(self.unallocatedRange.tolist(), len(plates), mode))

        jd0 = self.startDate

        while jd0 <= self.endDate:

            if useDateAtAPO:
                platesToSchedule = [plate for plate in plates
                                    if plate.dateAtAPO is None or
                                    plate.dateAtAPO <= jd0]
            else:
                platesToSchedule = plates

            optimalPlate, newExposures = logic.getOptimalPlate(
                platesToSchedule, [jd0, self.endDate], mode=mode, **kwargs)

            if optimalPlate is None:
                if jd0 < self.endDate - 2 * expTimeJD:
                    jd0 += expTimeJD
                    continue
                else:
                    break
            else:

                if optimalPlate in self.scheduled:
                    # Removes old version of optimalPlate in self.scheduled.
                    self.scheduled.remove(optimalPlate)
                self.scheduled.append(optimalPlate)

                self.allocateJDs(newExposures)
                jd0 = np.max([exp.getJD()[1] for exp in newExposures])

                # Defines some flags to be logged if plate is in Cosmic, being
                # drilled or not on the mountain.

                flags = ''

                if (optimalPlate.drilled is False and
                        not isinstance(optimalPlate, Field)):
                    flags = _color_text('** not yet drilled **', 'white')
                elif not isinstance(optimalPlate, Field):
                    location = optimalPlate.getLocation()  # May be None
                    if not location and not isinstance(optimalPlate, Field):
                        flags = _color_text('** unknown location **', 'red')
                    elif location != 'APO' and location != 'Cosmic':
                        flags = _color_text('** not on the mountain **',
                                            'red')
                    elif location == 'Cosmic':
                        flags = _color_text('** in Cosmic **', 'red')

                # Calculates completion before and after the simulation.
                completionPre = self.calculatePlateCompletion(
                    optimalPlate, rejectExposures=newExposures, useMock=True)
                completionPost = self.calculatePlateCompletion(
                    optimalPlate, useMock=True)
                nExps = len(newExposures)

                # Gets number of orphaned exposures after the simulation
                nOrphanedPost = 0
                for ss in optimalPlate.sets:
                    if ss.getStatus()[0] in ['Incomplete', 'Unplugged']:
                        for exp in ss.totoroExposures:
                            if exp.isMock:
                                nOrphanedPost += 1

                # Logs the result of the simulation. Format changed depending
                # on whether this is plate or a field.
                if not isinstance(optimalPlate, Field):
                    log.info('...... plate_id={0}, ({1} new exps, '
                             '{2:.2f} -> {3:.2f} complete) {4}'
                             .format(optimalPlate.plate_id,
                                     nExps, completionPre, completionPost,
                                     flags))

                    if nOrphanedPost > 0 and mode == 'plugger':
                        warnings.warn('... plate_id={0} has {1} orphaned '
                                      'exps after simulation'
                                      .format(optimalPlate.plate_id,
                                              nOrphanedPost),
                                      TotoroUserWarning)

                else:
                    log.info(_color_text(
                                '...... manga_tiledid={0} ({1} new exps, '
                                '{2:.2f} -> {3:.2f} complete) {4}'
                                .format(optimalPlate.getMangaTileID(),
                                        nExps, completionPre, completionPost,
                                        flags),
                                'yellow'))

        if showUnobservedTimes:
            if self.remainingTime <= expTimeJD:
                return True
            else:
                log.info('... unobserved times: {0}'.format(
                         str(self.unallocatedRange)))
                return False
