#!/usr/bin/env python
# encoding: utf-8
"""
planner.py

Created by José Sánchez-Gallego on 27 Oct 2014.
Licensed under a 3-clause BSD license.

Revision history:
    27 Oct 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.internal.manga.Totoro import log, config, readPath
from sdss.internal.manga.Totoro.dbclasses import Fields
from sdss.internal.manga.Totoro.scheduler.timeline import Timelines
from sdss.internal.manga.Totoro import exceptions
from astropy import time
import warnings


class PlannerScheduler(object):
    """A class for field selection."""

    def __init__(self, schedule, **kwargs):

        self.timelines = Timelines(schedule, **kwargs)
        log.debug('created PlannerScheduler with {0} timelines'
                  .format(len(self.timelines)))

        self.plates = self.getPlates(updateSets=False, silent=True, **kwargs)
        drilling = [plate for plate in self.plates if plate.isMock]
        log.info('Plates found: {0} ({1} in process of being drilled)'
                 .format(len(self.plates), len(drilling)))

        self.fields = Fields()

    def getPlates(self, **kwargs):
        """Gets plates that are already drilled, with some filtering."""

        from sdss.internal.manga.Totoro.dbclasses import (getAll, Plate,
                                                          getTilingCatalogue)

        plates = getAll(updateSets=False, silent=True)

        # Selects only non-started plates with priority > 1
        plates = [plate for plate in plates
                  if plate.getPlateCompletion() == 0.0 and plate.priority > 1
                  and len(plate.getTotoroExposures()) == 0]

        tilesBeingDrilled = map(
            int, open(
                readPath(config['planner']['tilesBeingDrilled']),
                'r').read().splitlines())

        tiles = getTilingCatalogue()

        for tile in tilesBeingDrilled:
            if tile not in tiles['ID']:
                warnings.warn('tile being drilled not identified, '
                              'manga_tileid={0}'.format(tile),
                              exceptions.TotoroUserWarning)
                continue
            row = tiles[tiles['ID'] == tile]
            plates.append(
                Plate.createMockPlate(ra=row['RA'][0], dec=row['DEC'][0],
                                      silent=True))

        platesIgnore = map(
            int,
            open(readPath(config['planner']['platesIgnore']),
                 'r').read().splitlines())
        plates = [plate for plate in plates
                  if plate.plate_id not in platesIgnore]

        return plates

    def schedule(self, **kwargs):
        """Runs the scheduling simulation."""

        from sdss.internal.manga.Totoro.dbclasses import Field

        for nn, timeline in enumerate(self.timelines):

            startDate = time.Time(timeline.startDate, format='jd')

            timeline.schedule(self.plates, mode='planner', **kwargs)

            remainingTime = timeline.remainingTime
            totalTime = 24. * (timeline.endDate - timeline.startDate)

            log.info('Scheduled timeline {0:.3f}-{1:.3f} [{2}] ({3:.1f}h). '
                     .format(timeline.startDate, timeline.endDate,
                             startDate.iso.split()[0], totalTime))
            log.info('... plates observed: {0} (Unused time {1:.1f}h)'
                     .format(len(timeline._plates), remainingTime))

            if remainingTime > 0:
                timeline.schedule(self.fields, mode='planner', **kwargs)
                nFields = [field for field in timeline._plates
                           if isinstance(field, Field)]
                remainingTime = timeline.remainingTime
                log.info('... fields observed: {0} (unused time {1:.1f}h)'
                         .format(len(nFields), remainingTime))
