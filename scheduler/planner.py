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
from sdss.internal.manga.Totoro.core.colourPrint import _color_text
from astropy import table
from astropy import time
import warnings
import numpy as np
import os


expTime = config['exposure']['exposureTime']
expTimeEff = expTime / config['planner']['efficiency']
minimumPlugPriority = config['planner']['noPlugPriority']


class PlannerScheduler(object):
    """A class for field selection."""

    def __init__(self, schedule, useFields=True, **kwargs):

        self.timelines = Timelines(schedule, **kwargs)
        log.debug('created PlannerScheduler with {0} timelines'
                  .format(len(self.timelines)))

        # Selects all valid plates, including complete ones.
        allPlates = self.getPlates(**kwargs)

        # Selects only plates that are incomplete and have enough priority
        self.plates = [plate for plate in allPlates
                       if len(plate.getTotoroExposures()) == 0 and
                       plate.priority > minimumPlugPriority]

        drilling = [plate for plate in self.plates if not plate.drilled]

        txtDrilling = ''
        if len(drilling) > 0:
            txtDrilling = _color_text(
                '({0} in process of being drilled)'.format(len(drilling)),
                'red')

        log.info('Plates found: {0} {1}'.format(len(self.plates), txtDrilling))

        # Gets fields (rejectDrilled=False because we do our own rejection)
        if useFields:
            self.createFields(allPlates)
        else:
            self.fields = []

    def createFields(self, allPlates):
        """Creates a field list from a tiling catalogue."""

        try:
            scienceCatalogue = table.Table.read(
                readPath(config['fields']['scienceCatalogue']))
            if 'MANGA_TILEID' not in scienceCatalogue.columns:
                warnings.warn('science catalogue does not contain '
                              'MANGA_TILEID. Won\'t check for number of '
                              'targets', exceptions.TotoroPlannerWarning)
                scienceCatalogue = None
        except:
            scienceCatalogue = None
            warnings.warn('science catalogue cannot be found. '
                          'Won\'t check for number of targets',
                          exceptions.TotoroPlannerWarning)

        from sdss.internal.manga.Totoro.dbclasses import Fields

        tmpFields = Fields(rejectDrilled=False)
        tileIDPlates = [plate.manga_tileid for plate in allPlates]

        self.fields = []

        for field in tmpFields:

            if field.manga_tileid in tileIDPlates:
                continue

            if scienceCatalogue is not None:

                scienceCatRows = scienceCatalogue[
                    scienceCatalogue['MANGA_TILEID'] == field.manga_tileid]

                if len(scienceCatRows) == 0:
                    log.debug('no targets for manga_tileid={0}. Skipping.'
                              .format(field.manga_tileid))
                    continue

            self.fields.append(field)

        nFieldsDrilled = len(tmpFields) - len(self.fields)
        if nFieldsDrilled > 0:
            log.info('rejected {0} fields because they have already '
                     'been drilled or have no targets.'.format(nFieldsDrilled))

    @staticmethod
    def getPlates(usePlatesNotAtAPO=True, **kwargs):
        """Gets plates that are already drilled or in process of being so,
        with some filtering."""

        from sdss.internal.manga.Totoro.dbclasses import (getAll, Plate,
                                                          getTilingCatalogue)

        allPlates = getAll(rejectSpecial=True, updateSets=False, silent=True,
                           fullCheck=False)

        # Selects plates with valid statuses
        validPlates = []
        for plate in allPlates:
            validPlate = True
            for status in plate.statuses:
                if status.label in ['Rejected', 'Unobservable']:
                    print(plate)
                    validPlate = False
            if (not usePlatesNotAtAPO and
                    plate.getLocation() not in ['APO', 'Cosmic']):
                validPlate = False
            if validPlate:
                validPlates.append(plate)

        # Adds tiles being drilled from file
        if ('tilesBeingDrilled' in config['fields'] and
                config['fields']['tilesBeingDrilled'].lower() != 'none' and
                usePlatesNotAtAPO):

            tilesPath = readPath(config['fields']['tilesBeingDrilled'])

            if not os.path.exists(tilesPath):

                warnings.warn(
                    'tilesBeingDrilled path does not exists. Skipping '
                    'additional tiles.', exceptions.TotoroPlannerWarning)

            else:

                tilesBeingDrilledRaw = open(tilesPath, 'r') \
                    .read().strip().splitlines()

                if (len(tilesBeingDrilledRaw) > 0
                        and tilesBeingDrilledRaw[0] != ''):
                    tilesBeingDrilled = map(int, tilesBeingDrilledRaw)
                else:
                    tilesBeingDrilled = []

                tiles = getTilingCatalogue()

                for tile in tilesBeingDrilled:

                    if tile not in tiles['ID']:
                        warnings.warn('tile being drilled not identified, '
                                      'manga_tileid={0}'.format(tile),
                                      exceptions.TotoroPlannerWarning)
                        continue

                    tileRow = tiles[tiles['ID'] == tile]

                    mockPlate = Plate.createMockPlate(
                        ra=tileRow['RA'][0], dec=tileRow['DEC'][0],
                        manga_tileid=tile, silent=True)

                    mockPlate.manga_tileid = tile
                    mockPlate.drilled = False

                    validPlates.append(mockPlate)

        # Plates to be rejected.
        if ('platesIgnore' in config['planner'] and
                config['planner']['platesIgnore'].lower() != 'none'):
            platesIgnore = map(
                int, open(readPath(config['planner']['platesIgnore']), 'r')
                .read().splitlines())
            validPlates = [plate for plate in validPlates
                           if plate.plate_id not in platesIgnore]

        return validPlates

    def schedule(self, useFields=True, goodWeatherFraction=None, **kwargs):
        """Runs the scheduling simulation."""

        goodWeatherFraction = goodWeatherFraction \
            if goodWeatherFraction is not None \
            else config['planner']['goodWeatherFraction']

        SN2_red = config['SN2thresholds']['plateRed']
        SN2_blue = config['SN2thresholds']['plateBlue']

        efficiency = kwargs.get('efficiency', config['planner']['efficiency'])

        log.info('Good weather fraction: {0:.2f}'.format(goodWeatherFraction))
        log.info('Efficiency: {0:.2f}'.format(efficiency))
        log.info('SN2 red={0:.1f}, blue={1:.1f}'.format(SN2_red, SN2_blue))

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
                log.info(
                    _color_text(
                        '... skipping timeline because of bad weather.',
                        'cyan'))
                timeline.observed = False
                continue

            timeline.observed = True

            timeline.schedule(self.plates, mode='planner', **kwargs)

            remainingTime = timeline.remainingTime
            log.info('... plates observed: {0} (Unused time {1:.1f}h)'
                     .format(len(timeline.plates), remainingTime))

            nAssigned += len(timeline.plates)

            if remainingTime > expTimeEff / 3600. and useFields:
                log.info(_color_text('now scheduling fields', 'yellow'))
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
                              exceptions.TotoroPlannerWarning)

    def getGoodWeatherIndices(self, goodWeatherFraction, seed=None, **kwargs):
        """Returns random indices with good weather."""

        np.random.seed(seed if seed is not None
                       else config['planner']['seed'])

        nTimelines = int(len(self.timelines) * goodWeatherFraction)
        indices = np.random.choice(np.arange(len(self.timelines)), nTimelines,
                                   replace=False)

        return np.sort(indices)
