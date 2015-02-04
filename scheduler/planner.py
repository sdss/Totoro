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
from sdss.internal.manga.Totoro.scheduler.timeline import Timelines
from sdss.internal.manga.Totoro import exceptions
from astropy import time
import warnings
import numpy as np


expTime = config['exposure']['exposureTime']
expTimeEff = expTime / config['planner']['efficiency']


class PlannerScheduler(object):
    """A class for field selection."""

    def __init__(self, schedule, useFields=True, **kwargs):

        from sdss.internal.manga.Totoro.dbclasses import Fields

        self.timelines = Timelines(schedule, **kwargs)
        log.debug('created PlannerScheduler with {0} timelines'
                  .format(len(self.timelines)))

        self.plates, self._completed = self.getPlates(
            updateSets=False, silent=True, **kwargs)

        drilling = [plate for plate in self.plates if plate.isMock]
        log.info('Plates found: {0} ({1} in process of being drilled)'
                 .format(len(self.plates), len(drilling)))

        # Gets fields (rejectDrilled=False because we do our own rejection)
        if useFields:
            tmpFields = Fields(rejectDrilled=False)
            tileIDPlates = [plate.getMangaTileID() for plate in self.plates]

            self.fields = [field for field in tmpFields
                           if field.manga_tileid not in tileIDPlates]
        else:
            self.fields = []

        nFieldsDrilled = len(tmpFields) - len(self.fields)
        if nFieldsDrilled > 0:
            log.info('rejected {0} fields because they have already '
                     'been drilled'.format(nFieldsDrilled))

    @staticmethod
    def getPlates(**kwargs):
        """Gets plates that are already drilled or in process of being so,
        with some filtering."""

        from sdss.internal.manga.Totoro.dbclasses import (getAll, Plate,
                                                          getTilingCatalogue)

        allPlates = getAll(updateSets=False, silent=True)

        # Selects only non-started plates with priority > minimum
        minimumPlugPriority = config['planner']['noPlugPriority']
        plates = [plate for plate in allPlates
                  if plate.getPlateCompletion() == 0.0
                  and plate.priority > minimumPlugPriority
                  and len(plate.getTotoroExposures()) == 0
                  and plate.plate_id > 7800]

        completed = [plate for plate in allPlates if plate.isComplete is True]

        # Adds tiles being drilled from file
        if ('tilesBeingDrilled' in config['planner'] and
                config['planner']['tilesBeingDrilled'].lower() != 'none'):

            tilesBeingDrilledRaw = open(
                readPath(config['planner']['tilesBeingDrilled']), 'r') \
                .read().splitlines()

            if len(tilesBeingDrilledRaw) > 0 and tilesBeingDrilledRaw[0] != '':
                tilesBeingDrilled = map(int, tilesBeingDrilledRaw)
            else:
                tilesBeingDrilled = []

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

        return plates, completed

    def schedule(self, useFields=True, goodWeatherFraction=None, **kwargs):
        """Runs the scheduling simulation."""

        goodWeatherFraction = goodWeatherFraction \
            if goodWeatherFraction is not None \
            else config['planner']['goodWeatherFraction']

        efficiency = kwargs.get('efficiency', config['planner']['efficiency'])

        log.info('Good weather fraction: {0:.2f}'.format(goodWeatherFraction))
        log.info('Efficiency: {0:.2f}'.format(efficiency))

        goodWeatherIdx = self.getGoodWeatherIndices(goodWeatherFraction,
                                                    **kwargs)

        for nn, timeline in enumerate(self.timelines):

            nAssigned = 0

            startDate = time.Time(timeline.startDate, format='jd')
            totalTime = 24. * (timeline.endDate - timeline.startDate)

            log.info('Scheduling timeline {0:.3f}-{1:.3f} [{2}] ({3:.1f}h). '
                     .format(timeline.startDate, timeline.endDate,
                             startDate.iso.split()[0], totalTime))

            if nn not in goodWeatherIdx:
                log.info('... skipping timeline because of bad weather.')
                timeline.observed = False
                continue

            timeline.observed = True

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

    def getGoodWeatherIndices(self, goodWeatherFraction, seed=None, **kwargs):
        """Returns random indices with good weather."""

        np.random.seed(seed if seed is not None
                       else config['planner']['seed'])

        nTimelines = int(len(self.timelines) * goodWeatherFraction)
        indices = np.random.choice(np.arange(len(self.timelines)), nTimelines,
                                   replace=False)

        return np.sort(indices)
