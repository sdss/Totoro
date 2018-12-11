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

from __future__ import division, print_function

import numpy as np

from Totoro import config, log, site, utils
from Totoro.core.colourPrint import _color_text
from Totoro.exceptions import TotoroUserWarning
from Totoro.scheduler import scheduler_utils as logic


expTimeJD = config['exposure']['exposureTime'] / 86400.


class Timelines(list):
    """A list of Timelines."""

    def __init__(self, observingBlocks, mode='planner'):

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
            self.unallocatedRange = utils.removeInterval(self.unallocatedRange, expJD, wrapAt=None)

        return

    def calculatePlateCompletion(self, plate, rejectExposures=[], useMock=True):
        """Calculates plate completion after rejecting exposures."""

        from Totoro.dbclasses import Plate, Set

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

        return np.sum(unallocatedRange[:, 1] - unallocatedRange[:, 0]) * 24.

    def schedule(self, plates, mode='plugger', useDateAtAPO=True, **kwargs):
        """Schedules a list of plates.

        Uses a list of `plates` to schedule all the available time in the
        timeline.

        Parameters
        ----------
        plates : list
            The list of `Totoro.Plate` instances to schedule.
        mode : str
            Either `'plugger'` or `'planner'`. Depending on the `mode`, the
            scheduling algorithm will make different choices for the list
            of optimal plates to be scheduled.
        useDateAtAPO : bool
            If True, only plates with `Totoro.Plate.dateAtAPO <= startDate`
            will be used. This assumes that the user has somehow added the
            `dateAtAPO` information to the `plates` before calling the method.
        kwargs : dict
            Additional parameters to be passed to `getOptimalPlate`.

        """

        mode = mode.lower()

        log.debug('scheduling LST range {0} using {1} plates, mode={2}'.format(
            self.unallocatedRange.tolist(), len(plates), mode))

        jd0 = self.startDate

        if useDateAtAPO:
            platesToSchedule = [
                plate for plate in plates
                if plate.dateAtAPO is None or plate.dateAtAPO <= self.startDate
            ]
        else:
            platesToSchedule = plates

        while jd0 <= self.endDate:

            # Finds the optimal plate for the range [jd0, self.endDate]
            optimalPlate, newExposures = logic.getOptimalPlate(
                platesToSchedule, [jd0, self.endDate], mode=mode, **kwargs)

            if optimalPlate is None:
                # If no optimal plate is found, moves jd0 by one exposure time.
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

                # Updates self.unallocatedRange with the new exposures.
                self.allocateJDs(newExposures)

                # Updates the current jd0
                jd0 = np.max([exp.getJD()[1] for exp in newExposures])

                # Logs the scheduled plate
                self.log(optimalPlate, newExposures, mode)

        if self.remainingTime <= expTimeJD:
            return True
        else:
            lst_range = site.localSiderealTime(self.unallocatedRange)
            log.info('... unobserved times: {0}'.format(str(self.unallocatedRange)))
            log.info('... unobserved LSTs: {!r}'.format(lst_range.tolist()))
            return False

    def log(self, plate, newExposures, mode):
        """Logs an schedule plate."""

        from Totoro.dbclasses import Field

        # Defines some flags to be logged if plate is in Cosmic, being
        # drilled or not on the mountain.

        flags = []

        if plate.drilled is False and not isinstance(plate, Field):
            flags.append(_color_text('not yet drilled', 'white'))
        elif not isinstance(plate, Field):
            location = plate.getLocation()
            if not location and not isinstance(plate, Field):
                flags.append(_color_text('unknown location', 'red'))
            elif location == 'Cosmic' or location == 'Storage':
                flags.append(_color_text('in {0}'.format(location), 'red'))
            elif location == 'APO':
                pass
            else:
                flags.append(_color_text('not on the mountain', 'red'))

        # Checks if the plate is a backup
        if (len(plate.statuses) > 0 and plate.statuses[0].label.lower() == 'backup'):
            flags.append(_color_text('backup plate', 'green'))

        if not plate.inFootprint:
            flags.append(_color_text('not in footprint', 'red'))

        flagsStr = '** {0} **'.format('; '.join(flags)) \
            if len(flags) > 0 else ''

        # Calculates completion before and after the simulation.
        completionPre = self.calculatePlateCompletion(
            plate, rejectExposures=newExposures, useMock=True)
        completionPost = self.calculatePlateCompletion(plate, useMock=True)
        nExps = len(newExposures)

        # Gets number of orphaned exposures after the simulation
        nOrphanedPost = 0
        for ss in plate.sets:
            if ss.getStatus()[0] == 'Incomplete':
                nOrphanedPost += len([exp for exp in ss.totoroExposures])

        # Logs the result of the simulation. Format changed depending
        # on whether this is plate or a field.
        if not isinstance(plate, Field):
            log.info('...... plate_id={0}, ({1} new exps, '
                     '{2:.2f} -> {3:.2f} complete) {4}'
                     .format(plate.plate_id, nExps, completionPre, completionPost, flagsStr))

            if nOrphanedPost > 0 and mode == 'plugger':
                log.warning(
                    '... plate_id={0} has {1} orphaned '
                    'exps after simulation'.format(plate.plate_id, nOrphanedPost),
                    TotoroUserWarning)

        else:
            log.info(
                _color_text(
                    '...... manga_tiledid={0} ({1} new exps, '
                    '{2:.2f} -> {3:.2f} complete) {4}'
                    .format(plate.getMangaTileID(), nExps, completionPre, completionPost,
                            flagsStr), 'yellow'))
