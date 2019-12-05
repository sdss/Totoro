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

from __future__ import division, print_function

import os

import numpy as np
from astropy import table, time
from astropy.io.ascii.core import InconsistentTableError

from Totoro import config, exceptions, log, readPath, site
from Totoro.core.colourPrint import _color_text
from Totoro.scheduler import observingPlan
from Totoro.scheduler.timeline import Timelines
from Totoro.utils import isPlateComplete
from Totoro.utils.intervals import getIntervalIntersectionLength


__all__ = ['Planner']

expTime = config['exposure']['exposureTime']
expTimeEff = expTime / config['planner']['efficiency']
minimumPlugPriority = config['planner']['noPlugPriority']

# Time allocated to IC342
_ic342_allocated = 0.0

# LST range for IC342. This range will be skipped.
IC342_range = config['specialPrograms']['IC342']['lst_range']


class Planner(object):
    """A class for field selection.

    This class is intended for multi-day simulations either for full survey
    layout simulations or for field selection for plate design purposes.

    Parameters
    ----------
    startDate : float or None
        The JD of the beginning of the simulation. If `None`, the current JD
        is used.
    endDate : float or None
        The JD of the end of the simulation. If `None`, the last JD in the
        official schedule will be used.
    useFields : bool
        If `True`, not drilled fields will be used when drilled plated are not
        enough to allocate all the available time.
    plates : list or None
        Either a list of `Totoro.Plate` to be used or `None`, in which case the
        list is determined by `Planner.getPlates`.
    fields : list or None
        Either a list of `Totoro.Field` to use or `None`, in which case the
        list of fields is determined by `Planner.getFields`.
    rejectBackup : bool
        If True, plates marked as Backup will be rejected and their fields will
        not be reselected.
    excludeStarted : bool
        If True, rejects plates that have been started (have valid exposures).
    kwargs : dict
        Additional arguments to be passed to `Planner.getPlates`.

    """

    def __init__(self,
                 startDate=None,
                 endDate=None,
                 useFields=True,
                 plates=None,
                 fields=None,
                 rejectBackup=False,
                 excludeStarted=False,
                 **kwargs):

        self.startDate = time.Time.now().jd if startDate is None else startDate
        self.endDate = (observingPlan.plan[-1]['JD1'] if endDate is None else endDate)
        self._useFields = useFields

        # Forces plates and fields to be both defined or None.
        assert (not plates and not fields) or (plates and fields), \
            'plates and fields must be both None or a list.'

        assert self.startDate < self.endDate, 'startDate > endDate'

        self.blocks = observingPlan.getObservingBlocks(self.startDate, self.endDate)

        assert len(self.blocks) >= 1, 'no observing blocks selected'

        self.timelines = Timelines(self.blocks)
        self.unallocatedJDs = []
        log.debug('PLANNER: created Planner instance with {0} timelines'.format(
            len(self.timelines)))

        # If plates is defined, we use that list of plates
        if plates is not None:
            self.plates = plates
        else:
            # Selects all valid plates, including complete ones.
            validPlates = self.getPlates(**kwargs)

            # Selects only plates that are incomplete and have enough priority
            # self.plates = [plate for plate in validPlates
            #                if plate.priority > minimumPlugPriority]
            self.plates = [plate for plate in validPlates
                           if not isPlateComplete(plate, write_apocomplete=False,
                                                  mark_complete=False)]

            if excludeStarted:
                self.plates = [
                    plate for plate in self.plates if len(plate.getTotoroExposures()) == 0
                ]

            # Outputs the number of plates.
            drilling = [plate for plate in self.plates if not plate.drilled]

            txtDrilling = ''
            if len(drilling) > 0:
                txtDrilling = _color_text(
                    '({0} in process of being drilled)'.format(len(drilling)), 'red')

            log.info('PLANNER: Plates found: {0} {1}'.format(len(self.plates), txtDrilling))

        # Gets fields
        if fields is not None:
            self.fields = fields
        elif useFields:
            self.getFields(self.plates)

        if rejectBackup:
            n_pre = len(self.plates)
            self.plates = [
                plate for plate in self.plates
                if 'Backup' not in [status.label for status in plate.statuses]
            ]
            log.info('PLANNER: rejected {} backup plates'.format(n_pre - len(self.plates)))

    def getFields(self, plates):
        """Creates a field list from the tiling catalogue.

        Returns the list of tiles (as `Totoro.Field` instances) to be used
        when no drilled plates are available to cover the scheduled time.

        Tiles with `manga_tileids` that have already been drilled are rejected,
        unless the plate was priority == 1. Additionally, if a tile contains
        fewer than `config.fields.minTargetsInTile` targets, we skip it.

        """

        from Totoro.dbclasses import Fields

        try:
            scienceCatalogue = table.Table.read(readPath(config['fields']['scienceCatalogue']))
            if 'MANGA_TILEID' not in scienceCatalogue.columns:
                log.warning(
                    'PLANNER: science catalogue does not contain '
                    'MANGA_TILEID. Will not check for number of '
                    'targets', exceptions.TotoroPlannerWarning)
                scienceCatalogue = None
        except Exception:
            scienceCatalogue = None
            log.warning(
                'PLANNER: science catalogue cannot be found. '
                'Will not check for number of targets', exceptions.TotoroPlannerWarning)

        tmpFields = Fields(rejectDrilled=True, acceptPriority1=True)
        tileIDPlates = [plate.manga_tileid for plate in plates]
        self.fields = []

        weights = table.Table.read(
            readPath('+data/tile_weight.dat'), format='ascii.commented_header')

        for field in tmpFields:

            # Even if we have already rejected drilled fields, we double check.
            # For instance, if we have added plates that are actually tiles
            # in the process of being drilled, we want to make sure we don't
            # use those tiles again here.
            if field.manga_tileid in tileIDPlates:
                continue

            field.ancillary_weight = 0.0 \
                if field.manga_tileid not in weights['manga_tileid'] else \
                weights[weights['manga_tileid'] ==
                        field.manga_tileid]['ancillary_weight'][0]

            if scienceCatalogue is not None:

                scienceCatRows = scienceCatalogue[scienceCatalogue['MANGA_TILEID'] ==
                                                  field.manga_tileid]

                if len(scienceCatRows) < config['fields']['minTargetsInTile']:
                    log.debug('PLANNER: no targets for manga_tileid={0}. Skipping.'
                              .format(field.manga_tileid))
                    continue

            self.fields.append(field)

        nFieldsNoTargets = len(tmpFields) - len(self.fields)
        if nFieldsNoTargets > 0:
            log.info('PLANNER: rejected {0} fields because '
                     'they have not enough targets.'.format(nFieldsNoTargets))

    @staticmethod
    def getPlates(usePlatesNotAtAPO=True, useTilesBeingDrilled=True):
        """Gets the list of plates to schedule.

        Returns a list of `Totoro.Plate` instances with plates that have a
        valid status. Plates marked as `Rejected` or `Unobservable` are
        rejected, as are those marked special for MaNGA (all-sky, star plates,
        etc.)

        If `config.dateAtAPO` is defined, the date at which the plate will be
        available at APO will be included. This information is used during
        scheduling to determine if a plate can be observed at a certain time.
        `config.dateAtAPO` must be the path to a plaintext file with the format
        `plate_id, dateAtAPO`.

        This is a staticmethod and can be called independently.

        Parameters
        ----------
        usePlatesNotAtAPO : bool
            If True, plates that are already drilled but not at APO (they may
            be at Cosmic or in transit) will be considered. This includes
            plates already in the DB but not yet drilled.
        useTilesBeingDrilled : bool
            If True, tiles that are in the process of being drilled but have
            not yet been added to the DB will be considered as drilled plates.
            The list of `manga_tileids` to be considered should be given in a
            file whose path is defined in `config.fields.tilesBeingDrilled`.
            The file must contain as many lines as tiles to be considered
            drilled, with the format `manga_tileid,  dateAtAPO`. `dateAtAPO`
            can be omitted (e.g., `6325,`), in which case the tile will be made
            available immediately.

        Returns
        -------
        result : list
            A list of `Totoro.Plate` that match the input requirements.

        """

        from Totoro.dbclasses import getAll

        # Gets a list with all the plates
        allPlates = getAll(rejectSpecial=True, updateSets=False, silent=True, fullCheck=False)

        # Selects plates with valid statuses
        validPlates = []
        for plate in allPlates:
            statuses = [status.label for status in plate.statuses]
            if 'Rejected' in statuses or 'Unobservable' in statuses:
                continue
            if not usePlatesNotAtAPO and plate.getLocation() != 'APO':
                continue
            validPlates.append(plate)

        # Adds information about when the plates will be at APO.
        if config['dateAtAPO'].lower() != 'none':

            if not os.path.exists(readPath(config['dateAtAPO'])):
                log.warning('PLANNER: dateAtAPO file does not exists.',
                            exceptions.TotoroPlannerWarning)

            else:
                try:
                    dateAtAPO = table.Table.read(
                        readPath(config['dateAtAPO']),
                        format='ascii',
                        delimiter=',',
                        names=['plateid', 'jd'])

                    for plate in validPlates:
                        row = dateAtAPO[dateAtAPO['plateid'] == plate.plate_id]
                        if len(row) > 0:
                            plate.dateAtAPO = row['jd'][0]
                        else:
                            plate.dateAtAPO = 0.

                except InconsistentTableError:
                    log.warning(
                        'PLANNER: dateAtAPO file exists but could not be read.'
                        ' It is probably empty.', exceptions.TotoroPlannerWarning)

        if useTilesBeingDrilled:
            platesBeingDrilled = Planner._getPlatesBeingDrilled()
        else:
            platesBeingDrilled = []

        return validPlates + platesBeingDrilled

    @staticmethod
    def _getPlatesBeingDrilled():
        """Returns a list of mock plates with the tiles being drilled."""

        from Totoro.dbclasses import Plate, getTilingCatalogue

        # Checks that the file exists and can be read
        if ('tilesBeingDrilled' not in config['fields'] or
                config['fields']['tilesBeingDrilled'].lower() == 'none'):
            return []

        path = readPath(config['fields']['tilesBeingDrilled'].lower())
        if not os.path.exists(path):
            log.warning('PLANNER: tilesBeingDrilled file does not exist.',
                        exceptions.TotoroPlannerWarning)
            return []

        try:
            tilesBeingDrilled = table.Table.read(path, format='ascii', delimiter=',')
        except InconsistentTableError:
            log.warning(
                'PLANNER: tilesBeingDrilled could not be read although it '
                'exists. Make sure the file is not empty.', exceptions.TotoroPlannerWarning)
            return []
        except Exception:
            raise exceptions.TotoroPlannerError(
                'PLANNER: unknown error while reading tilesBeingDrilled')

        tiles = getTilingCatalogue()

        platesBeingDrilled = []

        for manga_tileid, dateAtAPO in tilesBeingDrilled:

            if manga_tileid not in tiles['ID']:
                log.warning(
                    'PLANNER: manga_tileid={0}: tile being drilled '
                    'not in tiling catalogue.'.format(manga_tileid),
                    exceptions.TotoroPlannerWarning)
                continue

            tileRow = tiles[tiles['ID'] == manga_tileid]

            mockPlate = Plate.createMockPlate(
                ra=tileRow['RA'][0], dec=tileRow['DEC'][0], manga_tileid=manga_tileid, silent=True)

            mockPlate.manga_tileid = manga_tileid
            mockPlate.drilled = False

            if dateAtAPO is not np.ma.masked:
                mockPlate.dateAtAPO = dateAtAPO

            platesBeingDrilled.append(mockPlate)

        return platesBeingDrilled

    def schedule(self,
                 goodWeatherFraction=config['planner']['goodWeatherFraction'],
                 efficiency=config['planner']['efficiency'],
                 prioritiseAPO=False,
                 useNightlyPlugged=True,
                 **kwargs):
        """Runs the scheduling simulation.

        Parameters
        ----------
        goodWeatherFraction : float
            The fraction of good weather to use. Defaults to
            `config.planner.goodWeatherFraction`.
        efficiency : float
            The efficiency to use to account for the overheads.Defaults to
            `config.planner.efficiency`.
        prioritiseAPO : bool
            If True, tries to find valid plates at APO first.
        useNightlyPlugged : bool
            If `True`, prioritises plates observed during the night (i.e.,
            simulates the concept of "plugged" plates).
        kwargs : dict
            Additional arguments to be passed to `getOptimalPlate`.

        """

        global _ic342_allocated

        goodWeatherFraction = goodWeatherFraction \
            if goodWeatherFraction is not None \
            else config['planner']['goodWeatherFraction']

        SN2_red = config['SN2thresholds']['plateRed']
        SN2_blue = config['SN2thresholds']['plateBlue']

        log.info('PLANNER: Good weather fraction: {0:.2f}'.format(goodWeatherFraction))
        log.info('PLANNER: Efficiency: {0:.2f}'.format(efficiency))
        log.info('PLANNER: SN2 red={0:.1f}, blue={1:.1f}'.format(SN2_red, SN2_blue))
        log.info('PLANNER: prioritise APO={0}'.format(prioritiseAPO))

        # Gets the indices of the timelines with good weather.
        goodWeatherIdx = self.getGoodWeatherIndices(goodWeatherFraction)

        for nn, timeline in enumerate(self.timelines):

            startDate = time.Time(timeline.startDate, format='jd')
            totalTime = 24. * (timeline.endDate - timeline.startDate)

            log.info('Scheduling timeline '
                     '{0:.3f}-{1:.3f} ({2:.2f}-{3:.2f}) [{4}] ({5:.1f}h). '.format(
                         timeline.startDate, timeline.endDate,
                         site.localSiderealTime(timeline.startDate),
                         site.localSiderealTime(timeline.endDate),
                         startDate.iso.split()[0], totalTime))

            ic342_lenght = getIntervalIntersectionLength(
                [site.localSiderealTime(timeline.startDate),
                 site.localSiderealTime(timeline.endDate)], IC342_range)

            _ic342_allocated += ic342_lenght

            if nn not in goodWeatherIdx:
                log.info(_color_text('... skipping timeline because of ' 'bad weather.', 'cyan'))
                timeline.observed = False
                continue

            timeline.observed = True

            if not self._useFields:
                timeline.schedule(
                    self.plates,
                    mode='planner',
                    prioritiseAPO=prioritiseAPO,
                    useNightlyPlugged=useNightlyPlugged,
                    **kwargs)
            else:
                timeline.schedule(
                    self.plates + self.fields,
                    mode='planner',
                    prioritiseAPO=prioritiseAPO,
                    useNightlyPlugged=useNightlyPlugged,
                    **kwargs)

            remainingTime = timeline.remainingTime
            if remainingTime < 0:
                remainingTime = 0.0

            colour = 'red' if remainingTime > 0.1 else 'default'
            log.info(
                _color_text(
                    '... plates observed: {0} (Unused time {1:.2f}h)'.format(
                        len(timeline.scheduled), remainingTime), colour))

            if remainingTime > 0:
                unallocatedThisTimeline = np.atleast_2d(timeline.unallocatedRange)
                for j0, j1 in unallocatedThisTimeline:
                    self.unallocatedJDs.append([j0, j1])

            nCarts = len(config['mangaCarts']) - len(config['offlineCarts'])
            if len(timeline.scheduled) > nCarts:
                log.warning(
                    'more plates ({0}) scheduled than carts available ({1})'.format(
                        len(timeline.scheduled), nCarts), exceptions.TotoroPlannerWarning)

        self.unallocatedJDs = np.array(self.unallocatedJDs)

        log.warning('IC342 allocated time {0:.1f}h'.format(_ic342_allocated),
                    exceptions.TotoroPlannerWarning)

    def getGoodWeatherIndices(self, goodWeatherFraction, seed=None):
        """Returns random indices with good weather."""

        np.random.seed(seed if seed is not None else config['planner']['seed'])

        nTimelines = int(len(self.timelines) * goodWeatherFraction)
        indices = np.random.choice(np.arange(len(self.timelines)), nTimelines, replace=False)

        return np.sort(indices)
