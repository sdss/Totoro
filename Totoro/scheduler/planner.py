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
from Totoro import log, config, readPath, site
from Totoro.scheduler.timeline import Timelines
from Totoro.scheduler import observingPlan
from Totoro.core.colourPrint import _color_text
from Totoro import exceptions
from astropy import table
from astropy import time
import warnings
import numpy as np
import os


__all__ = ['Planner']

expTime = config['exposure']['exposureTime']
expTimeEff = expTime / config['planner']['efficiency']
minimumPlugPriority = config['planner']['noPlugPriority']


class Planner(object):
    """A class for field selection."""

    def __init__(self, startDate=None, endDate=None, useFields=True,
                 optimiseFootprint=True, plates=None, fields=None, **kwargs):

        log.info('entering PLANNER mode.')

        if startDate is None:
            startDate = time.Time.now().jd
        if endDate is None:
            endDate = observingPlan.plan[-1]['JD1']

        self.startDate = startDate
        self.endDate = endDate

        assert self.startDate < self.endDate, 'startDate > endDate'

        self.blocks = observingPlan.getObservingBlocks(self.startDate,
                                                       self.endDate)

        assert len(self.blocks) >= 1, 'no observing blocks selected'

        self.timelines = Timelines(self.blocks, **kwargs)
        log.debug('created PlannerScheduler with {0} timelines'
                  .format(len(self.timelines)))

        # If plates is defined, we use that list of plates
        if plates is not None:
            self.plates = plates
        else:
            # Selects all valid plates, including complete ones.
            self._allPlates = self.getPlates(
                optimiseFootprint=optimiseFootprint, **kwargs)

            # Selects only plates that are incomplete and have enough priority
            self.plates = [plate for plate in self._allPlates
                           if len(plate.getTotoroExposures()) == 0 and
                           plate.priority > minimumPlugPriority]

            drilling = [plate for plate in self.plates if not plate.drilled]

            txtDrilling = ''
            if len(drilling) > 0:
                txtDrilling = _color_text(
                    '({0} in process of being drilled)'.format(len(drilling)),
                    'red')

            log.info('Plates found: {0} {1}'.format(len(self.plates),
                                                    txtDrilling))

        # Gets fields (rejectDrilled=False because we do our own rejection)
        self.fields = []
        if fields is not None:
            self.fields = fields
        elif useFields:
            self.createFields(self._allPlates)

        # Defines priorities
        if optimiseFootprint:
            self._assignPriorities()

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

        from Totoro.dbclasses import Fields

        tmpFields = Fields(rejectDrilled=False)
        tileIDPlates = [plate.manga_tileid for plate in allPlates]

        self.fields = []

        for field in tmpFields:

            if field.manga_tileid in tileIDPlates:
                continue

            if scienceCatalogue is not None:

                scienceCatRows = scienceCatalogue[
                    scienceCatalogue['MANGA_TILEID'] == field.manga_tileid]

                if len(scienceCatRows) < 12:
                    log.debug('no targets for manga_tileid={0}. Skipping.'
                              .format(field.manga_tileid))
                    continue

            self.fields.append(field)

        nFieldsDrilled = len(tmpFields) - len(self.fields)
        if nFieldsDrilled > 0:
            log.info('rejected {0} fields because they have already '
                     'been drilled or have no targets.'.format(nFieldsDrilled))

    @staticmethod
    def getPlates(usePlatesNotAtAPO=True, usePlatesBeingDrilled=True,
                  optimiseFootprint=True, skipDrilled=False, **kwargs):
        """Gets plates that are already drilled or in process of being so,
        with some filtering."""

        from Totoro.dbclasses import getAll, Plate, getTilingCatalogue

        allPlates = getAll(rejectSpecial=True, updateSets=False, silent=True,
                           fullCheck=False)

        # Selects plates with valid statuses
        validPlates = []
        for plate in allPlates:
            validPlate = True
            for status in plate.statuses:
                if status.label in ['Rejected', 'Unobservable']:
                    validPlate = False
            if not usePlatesNotAtAPO and plate.getLocation() != 'APO':
                validPlate = False
            if optimiseFootprint:
                if ((plate.ra < 100 or plate.ra > 300) and
                        (plate.dec < -1 or plate.dec > 1)):
                    validPlate = False
            if skipDrilled:
                if len(plate.getScienceExposures()) == 0:
                    validPlate = False
            if validPlate:
                validPlates.append(plate)

        # Now we open the dateAtAPO file and add the date to the corresponging
        # plates. This is useful when, even if the plate is at APO, you want to
        # delay its scheduling until a certain time.
        if config['dateAtAPO'].lower() != 'none':

            dateAtAPO_Path = readPath(config['dateAtAPO'])

            if not os.path.exists(dateAtAPO_Path):
                warnings.warn('dateAtAPO file does not exists.',
                              exceptions.TotoroPlannerWarning)

            else:
                dateAtAPO_Table = table.Table.read(
                    dateAtAPO_Path, format='ascii.commented_header',
                    delimiter=',')

                for plate_id, dateAtAPO in dateAtAPO_Table:
                    for plate in validPlates:
                        if plate.plate_id == plate_id:
                            plate.dateAtAPO = dateAtAPO

        # Adds tiles being drilled from file
        if (config['fields']['tilesBeingDrilled'].lower() != 'none' and
                usePlatesBeingDrilled):

            tilesPath = readPath(config['fields']['tilesBeingDrilled'])

            if not os.path.exists(tilesPath):
                warnings.warn(
                    'tilesBeingDrilled path does not exists. Skipping '
                    'additional tiles.', exceptions.TotoroPlannerWarning)

            else:

                tilesBeingDrilled = table.Table.read(
                    tilesPath, format='ascii.commented_header', delimiter=',')

                if len(tilesBeingDrilled) > 0:

                    tiles = getTilingCatalogue()

                    for manga_tileid, dateAtAPO in tilesBeingDrilled:

                        if manga_tileid not in tiles['ID']:
                            warnings.warn('manga_tileid={0}: tile being '
                                          'drilled not in tiling '
                                          'catalogue.'.format(manga_tileid),
                                          exceptions.TotoroPlannerWarning)
                            continue

                        tileRow = tiles[tiles['ID'] == manga_tileid]

                        mockPlate = Plate.createMockPlate(
                            ra=tileRow['RA'][0], dec=tileRow['DEC'][0],
                            manga_tileid=manga_tileid, silent=True)

                        mockPlate.manga_tileid = manga_tileid
                        mockPlate.drilled = False

                        if dateAtAPO is not np.ma.masked:
                            mockPlate.dateAtAPO = dateAtAPO

                        validPlates.append(mockPlate)

        return validPlates

    def _assignPriorities(self):
        """Defines footprint priorities."""

        for field in self.fields:
            if ((field.ra < 100 or field.ra > 300) and
                    (field.dec > -1 and field.dec < 1)):
                field.priority = 9
            elif ((field.ra > 199.5 - 7.5 and field.ra < 199.5 + 7.5) and
                    (field.dec > 29 - 5 and field.dec < 29 + 5)):
                field.priority = 9
            elif ((field.ra > 11 * 15. and field.ra < 14 * 15.) and
                    (field.dec > 45 and field.dec < 50)):
                field.priority = 9
            elif field.dec < 20:
                field.priority = 6.5

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

            log.info('Scheduling timeline {0:.3f}-{1:.3f} ({2:.2f}-{3:.2f}) '
                     '[{4}] ({5:.1f}h). '
                     .format(timeline.startDate, timeline.endDate,
                             site.localSiderealTime(timeline.startDate),
                             site.localSiderealTime(timeline.endDate),
                             startDate.iso.split()[0], totalTime))

            if nn not in goodWeatherIdx:
                log.info(_color_text('... skipping timeline because of '
                                     'bad weather.', 'cyan'))
                timeline.observed = False
                continue

            timeline.observed = True

            timeline.schedule(self.plates + self.fields, mode='planner',
                              **kwargs)

            remainingTime = timeline.remainingTime
            colour = 'red' if remainingTime > 0.1 else 'default'
            log.info(
                _color_text(
                    '... plates observed: {0} (Unused time {1:.2f}h)'
                    .format(len(timeline.scheduled), remainingTime), colour))

            nAssigned += len(timeline.scheduled)

            nCarts = len(config['mangaCarts']) - len(config['offlineCarts'])
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
