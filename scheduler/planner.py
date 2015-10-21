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


expTime = config['exposure']['exposureTime']
expTimeEff = expTime / config['planner']['efficiency']


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

        # Gets fields (rejectDrilled=False because we do our own rejection)
        tmpFields = Fields(rejectDrilled=False)
        tileIDPlates = [plate.getMangaTileID() for plate in self.plates]

        self.fields = [field for field in tmpFields
                       if field.manga_tileid not in tileIDPlates]

        nFieldsDrilled = len(tmpFields) - len(self.fields)
        if nFieldsDrilled > 0:
            log.info('rejected {0} fields because they have already '
                     'been drilled'.format(nFieldsDrilled))

    def getPlates(self, **kwargs):
        """Gets plates that are already drilled or in process of being so,
        with some filtering."""

        from sdss.internal.manga.Totoro.dbclasses import (getAll, Plate,
                                                          getTilingCatalogue)

        plates = getAll(updateSets=False, silent=True)

        # Selects only non-started plates with priority > minimum
        minimumPlugPriority = config['planner']['noPlugPriority']
        plates = [plate for plate in plates
                  if plate.getPlateCompletion() == 0.0
                  and plate.priority > minimumPlugPriority
                  and len(plate.getTotoroExposures()) == 0
                  and plate.plate_id > 7800]

        # Adds tiles being drilled from file
        if ('tilesBeingDrilled' in config['planner'] and
                config['planner']['tilesBeingDrilled'].lower() != 'none'):

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
                mockPlate = Plate.createMockPlate(
                    ra=row['RA'][0], dec=row['DEC'][0], silent=True)
                mockPlate.manga_tileid = tile
                plates.append(mockPlate)

        # Plates to be rejected.
        if ('platesIgnore' in config['planner'] and
                config['planner']['platesIgnore'].lower() != 'none'):
            platesIgnore = map(
                int, open(readPath(config['planner']['platesIgnore']), 'r')
                .read().splitlines())
            plates = [plate for plate in plates
                      if plate.plate_id not in platesIgnore]

        return plates

    def schedule(self, useFields=True, **kwargs):
        """Runs the scheduling simulation."""

        for nn, timeline in enumerate(self.timelines):

            nAssigned = 0

            startDate = time.Time(timeline.startDate, format='jd')
            totalTime = 24. * (timeline.endDate - timeline.startDate)

            log.info('Scheduling timeline {0:.3f}-{1:.3f} [{2}] ({3:.1f}h). '
                     .format(timeline.startDate, timeline.endDate,
                             startDate.iso.split()[0], totalTime))

            timeline.schedule(self.plates, mode='planner', **kwargs)

            remainingTime = timeline.remainingTime
            log.info('... plates observed: {0} (Unused time {1:.1f}h)'
                     .format(len(timeline.plates), remainingTime))

            nAssigned += len(timeline.plates)

            if remainingTime > expTimeEff / 3600. and useFields:
                log.info('now scheduling fields')
                timeline.schedule(self.fields, mode='planner', **kwargs)
                nFields = len(timeline.plates) - nAssigned
                remainingTime = timeline.remainingTime
                log.info('... fields observed: {0} (unused time {1:.1f}h)'
                         .format(nFields, remainingTime))
                nAssigned += nFields

            nCarts = len(config['carts'])
            if nAssigned > nCarts:
                warnings.warn('more plates ({0}) scheduled than carts '
                              'available ({1})'.format(nAssigned, nCarts),
                              exceptions.TotoroWarning)
