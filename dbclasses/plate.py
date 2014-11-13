#!/usr/bin/env python
# encoding: utf-8
"""
plate.py

Created by José Sánchez-Gallego on 18 Apr 2014.
Licensed under a 3-clause BSD license.

Revision history:
    18 Apr 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.internal.manga.Totoro.apoDB import TotoroDBConnection
from sdss.internal.manga.Totoro import utils
from sdss.internal.manga.Totoro import exceptions as TotoroExpections
from sdss.internal.manga.Totoro import logic
from sdss.internal.manga.Totoro import log, config, dustMap, site
# from sdss.internal.manga.Totoro import scheduler
import warnings
from astropy import time
import set as TotoroSet
from exposure import Exposure
import numpy as np
from copy import deepcopy
import os


__ALL__ = ['getPlugged', 'getAtAPO', 'getAll', 'Plates', 'Plate',
           'fromPlateID']


totoroDB = TotoroDBConnection()
plateDB = totoroDB.plateDB
mangaDB = totoroDB.mangaDB


def getPlugged(onlyIncomplete=False, **kwargs):

    session = totoroDB.Session()
    with session.begin(subtransactions=True):
        activePluggings = session.query(
            plateDB.ActivePlugging).join(
                plateDB.Plugging,
                plateDB.Plate,
                plateDB.PlateToSurvey,
                plateDB.Survey,
                plateDB.SurveyMode).filter(
                    plateDB.Survey.label == 'MaNGA',
                    plateDB.SurveyMode.label == 'MaNGA dither').order_by(
                        plateDB.Plate.plate_id).all()

    plates = [actPlug.plugging.plate_pk for actPlug in activePluggings]

    if onlyIncomplete:
        return _getIncomplete(plates, **kwargs)
    else:
        return Plates(plates, **kwargs)


def getAtAPO(onlyIncomplete=False, **kwargs):

    session = totoroDB.Session()
    with session.begin(subtransactions=True):
        plates = session.query(plateDB.Plate).join(
            plateDB.PlateToSurvey, plateDB.Survey, plateDB.SurveyMode).filter(
                plateDB.Survey.label == 'MaNGA',
                plateDB.SurveyMode.label == 'MaNGA dither'
            ).join(plateDB.PlateLocation).filter(
                plateDB.PlateLocation.label == 'APO').order_by(
                    plateDB.Plate.plate_id).all()

    plates = [plate.pk for plate in plates]

    if onlyIncomplete:
        return _getIncomplete(plates, **kwargs)
    else:
        return Plates(plates, **kwargs)


def getAll(onlyIncomplete=False, **kwargs):

    session = totoroDB.Session()
    with session.begin(subtransactions=True):
        plates = session.query(plateDB.Plate).join(
            plateDB.PlateToSurvey, plateDB.Survey, plateDB.SurveyMode
            ).filter(plateDB.Survey.label == 'MaNGA',
                     plateDB.SurveyMode.label == 'MaNGA dither').order_by(
                plateDB.Plate.plate_id).all()

    plates = [plate.pk for plate in plates]

    if onlyIncomplete:
        return _getIncomplete(plates, **kwargs)
    else:
        return Plates(plates, **kwargs)


def _getIncomplete(plates, **kwargs):

    totoroPlates = [Plate(plate, format='pk') for plate in plates]

    incompletePlates = []
    for totoroPlate in totoroPlates:
        if not utils.isPlateComplete(totoroPlate):
            incompletePlates.append(totoroPlate)
    return Plates(incompletePlates)


class Plates(list):

    def __init__(self, inp, format='pk', **kwargs):

        if all([isinstance(ii, Plate) for ii in inp]):
            list.__init__(self, inp)
        else:
            list.__init__(self, [Plate(ii, format=format, **kwargs)
                                 for ii in inp])

    @staticmethod
    def getPlugged(**kwargs):
        """For backards compatibility."""
        return getPlugged(**kwargs)

    @staticmethod
    def getAll(**kwargs):
        """For backards compatibility."""
        return getPlugged(**kwargs)

    @staticmethod
    def getAtAPO(**kwargs):
        """For backards compatibility."""
        return getPlugged(**kwargs)


def fromPlateID(plateid, **kwargs):
    return Plate(plateid, format='plate_id', **kwargs)


class Plate(plateDB.Plate):

    def __new__(cls, input=None, format='pk', **kwargs):

        if input is None:
            return plateDB.Plate.__new__(cls)

        base = cls.__bases__[0]

        if isinstance(input, base):
            instance = input
        else:
            session = totoroDB.Session()
            with session.begin(subtransactions=True):
                instance = session.query(base).filter(
                    eval('{0}.{1} == {2}'.format(base.__name__, format, input))
                    ).one()

        instance.__class__ = cls

        return instance

    def __init__(self, input, format='pk', mock=False, silent=False,
                 updateSets=True, mjd=None, **kwargs):

        self._complete = None
        self.isMock = mock
        self._kwargs = kwargs
        self.mjd = mjd

        if 'dust' in kwargs:
            self.dust = kwargs['dust']
        else:
            self.dust = None if dustMap is None else dustMap(self.ra, self.dec)

        self.mlhalimit = utils.mlhalimit(self.dec)

        if not self.isMock:
            self.checkPlate()
            self.sets = TotoroSet.getPlateSets(
                self.pk, format='pk', silent=silent, **kwargs)

            if silent is False:
                log.debug('loaded plate with pk={0}, plateid={1}'.format(
                          self.pk, self.plate_id))

            if updateSets:
                self.updatePlate()

        else:
            self.sets = []

    def __repr__(self):
        return ('<Totoro Plate (plate_id={0}, pk={1}, completion={2:.2f})>'
                .format(self.plate_id, self.pk, self.getPlateCompletion()))

    @classmethod
    def fromSets(cls, sets, silent=False, **kwargs):

        newPlate = plateDB.Plate.__new__(cls)
        newPlate.__init__(None, mock=True, **kwargs)
        newPlate.sets = sets

        if not silent:
            log.debug('created mock plate from sets pk={0}'.format(
                ', '.join(map(str, sets))))

        return newPlate

    @classmethod
    def createMockPlate(cls, ra=None, dec=None, silent=False, **kwargs):

        if ra is None or dec is None:
            raise TotoroExpections.TotoroError('ra and dec must be specified')

        mockPlate = plateDB.Plate.__new__(cls)
        mockPlate.__init__(None, mock=True, ra=ra, dec=dec, **kwargs)

        mockPlate.isMock = True

        if not silent:
            log.debug(
                'created mock plate with ra={0:.3f} and dec={0:.3f}'.format(
                    mockPlate.coords[0], mockPlate.coords[1]))

        return mockPlate

    def addExposure(self, exposure):
        return logic.addExposure(exposure, self)

    def checkPlate(self):

        if not hasattr(self, 'pk') or not hasattr(self, 'plate_id'):
            raise AttributeError('Plate instance has no pk or plate_id.')

        if not self.isMaNGA:
            raise TotoroExpections.TotoroError('this is not a MaNGA plate!')

        nMaNGAExposures = len(self.getMangadbExposures())
        nScienceExposures = len(self.getScienceExposures())
        if nMaNGAExposures != nScienceExposures:
            warnings.warn('plate_id={1}: {0} plateDB.Exposures found '
                          'but only {2} mangaDB.Exposures'.format(
                              nScienceExposures, self.plate_id,
                              nMaNGAExposures),
                          TotoroExpections.NoMangaExposure)

    @property
    def isMaNGA(self):

        for survey in self.surveys:
            if survey.label == 'MaNGA':
                return True

        return False

    def updatePlate(self, force=False):

        if self.isComplete:
            if not force:
                log.debug('plate_id={0} is marked complete. Not updating sets.'
                          .format(self.plate_id))
                return False
            else:
                log.debug('plate_id={0} is marked complete but force=True.'
                          .format(self.plate_id))

        if (self.currentSurveyMode is None or
                self.currentSurveyMode.label != 'MaNGA dither'):
            log.debug('plate_id={0} has surveyMode which is not MaNGA dither. '
                      'Not updating sets.')
            return False

        result = logic.updatePlate(self)
        if result:
            self.update()

        return True

    def rearrangeSets(self, startDate=None, **kwargs):
        result = logic.setArrangement.rearrangeSets(
            self, startDate=startDate, **kwargs)
        if result is True:
            self.update()

    @property
    def ra(self):
        return self.getCoordinates()[0]

    @property
    def dec(self):
        return self.getCoordinates()[1]

    @property
    def coords(self):
        """Alias for getCoordinates(). Provides backwards compatibility."""
        return self.getCoordinates()

    def getCoordinates(self):

        if 'ra' in self._kwargs and 'dec' in self._kwargs:
            if (self._kwargs['ra'] is not None and
                    self._kwargs['dec'] is not None):
                return np.array(
                    [self._kwargs['ra'], self._kwargs['dec']], np.float)
        else:

            if len(self.plate_pointings) > 1:
                warnings.warn('plate_id={0:d}: multiple plate pointings found.'
                              ' Using the first one.'.format(self.plate_id),
                              TotoroExpections.MultiplePlatePointings)

            return np.array(
                [self.plate_pointings[0].pointing.center_ra,
                 self.plate_pointings[0].pointing.center_dec], np.float)

    @property
    def isComplete(self):
        if self._complete is not None:
            return self._complete
        else:
            return utils.isPlateComplete(self)

    def copy(self):
        return deepcopy(self)

    def getPlateCompletion(self, includeIncompleteSets=False):

        totalSN = self.getCumulatedSN2(
            includeIncomplete=includeIncompleteSets)

        completionRates = totalSN.copy()
        completionRates[0:2] /= config['SN2thresholds']['plateBlue']
        completionRates[2:] /= config['SN2thresholds']['plateRed']

        return np.min(completionRates)

    def getCumulatedSN2(self, includeIncomplete=False):

        validStatuses = ['Good', 'Excellent']
        if includeIncomplete:
            validStatuses.append('Incomplete')

        validSets = []
        for set in self.sets:
            if set.getQuality()[0] in validStatuses:
                validSets.append(set)

        if len(validSets) == 0:
            return np.array([0.0, 0.0, 0.0, 0.0])
        else:
            return np.sum([set.getSN2Array() for set in validSets], axis=0)

    def getActiveCartNumber(self):

        for plugging in self.pluggings:
            if len(plugging.activePlugging) > 0:
                return int(plugging.cartridge.number)

        raise TotoroExpections.PlateNotPlugged(
            'plate_id={0} is not currently plugged'.format(self.plate_id))

    def getMangadbExposures(self):
        """Returns a list of mangaDB.Exposure objects with all the exposures
        for this plate."""

        session = totoroDB.Session()
        with session.begin(subtransactions=True):
            mangaExposures = session.query(totoroDB.mangaDB.Exposure).join(
                plateDB.Exposure,
                plateDB.Observation,
                plateDB.Plugging,
                plateDB.Plate).filter(
                    plateDB.Plate.pk == self.pk).all()

        return mangaExposures

    def getScienceExposures(self):
        """Returns a list of all plateDB.Exposure science exposures for
        the plate."""

        scienceExps = []

        for plugging in self.pluggings:
            scienceExps += plugging.scienceExposures()

        return scienceExps

    def getTotoroExposures(self):
        """Returns a list of Totoro dbclasses.Exposure instances for this
        Plate insance."""

        totoroExposures = []
        for set in self.sets:
            totoroExposures += set.totoroExposures

        return totoroExposures

    def getValidSets(self, includeIncomplete=False):

        validSets = []
        for set in self.sets:
            quality = set.getQuality()
            if quality in ['Good', 'Excellent', 'Poor']:
                validSets.append(set)
            elif includeIncomplete and quality == 'Incomplete':
                validSets.append(set)

        return validSets

    def update(self, **kwargs):
        """Reloads the plate."""

        kwargs.update(self._kwargs)

        newSelf = Plate(self.pk, format='pk', silent=True, mjd=self.mjd,
                        **kwargs)
        self = newSelf

        log.debug('plate_id={0} has been reloaded'.format(self.plate_id))

    def getValidExposures(self):
        """Returns all valid exposures, even if they belong to an incomplete
        or bad set."""

        validExposures = []

        for set in self.sets:
            for exp in set.totoroExposures:
                if exp.valid is True:
                    validExposures.append(exp)

        return validExposures

    def getHARange(self, intersect=False, mjd=None, **kwargs):

        ha0 = -self.mlhalimit
        ha1 = self.mlhalimit

        haRange = np.array([ha0, ha1]) % 360.
        haRange[haRange > 180.] -= 360.

        # if intersect is False:
        return haRange

        # if mjd is None:
        #     mjd = int(np.round(time.Time.now().mjd))
        #     log.debug('Plate.getHARange: MJD set to {0:d}'.format(mjd))
        # else:
        #     mjd = int(mjd)

        # jdRange = scheduler.observingPlan.getMJD(mjd)

        # if jdRange is None:
        #     warnings.warn('no observing block found for MJD={0:d}. '
        #                   'Observing windows will not be contrained.'
        #                   .format(mjd), TotoroExpections.NoObservingBlock)
        #     return haRange

        # observingRangeLST = np.array(map(site.localSiderealTime, jdRange))
        # observingRangeHA = (observingRangeLST * 15 - self.ra) % 360.

        # return utils.getIntervalIntersection(haRange, observingRangeHA)

    def getLSTRange(self, **kwargs):

        haRange = self.getHARange(**kwargs)

        ha0, ha1 = haRange

        lst0 = (ha0 + self.coords[0]) % 360. / 15
        lst1 = (ha1 + self.coords[0]) % 360. / 15

        LSTRange = np.array([lst0, lst1])

        return LSTRange

    def getUTVisibilityWindow(self, **kwargs):
        return self.getUTRange(**kwargs)

    def getUTRange(self, mjd=None, returnType='str', **kwargs):

        lst0, lst1 = self.getLSTRange(**kwargs)

        if mjd is None:
            mjd = int(np.round(time.Time.now().mjd))

        date = time.Time(mjd, format='mjd', scale='tai')
        date0, date1 = site.localSiderealTimeToDate([lst0, lst1], date=date)

        # Sanity check
        if date0.jd > date1.jd:
            raise ValueError('date0 is greater than date1')

        if returnType == 'str':
            return ('{0:%H:%M}'.format(date0.datetime),
                    '{0:%H:%M}'.format(date1.datetime))
        elif returnType == 'datetime':
            return (date0.datetime, date1.datetime)
        else:
            return (date0, date1)

    def isVisible(self, LST0, LST1, minLength=None):
        """Returns True if the plate is visible in a range of LST."""

        if minLength is None:
            minLength = config['exposure']['exposureTime']

        lstIntersectionLength = utils.getIntervalIntersectionLength(
            self.getLSTRange(), (LST0, LST1), wrapAt=24.)

        secIntersectionLength = lstIntersectionLength * 3600.

        return True if secIntersectionLength >= minLength else False

    def addMockExposure(self, exposure=None, startTime=None, set=None,
                        expTime=None, silent=False, **kwargs):
        """Creates a mock expusure in the best possible way."""

        ra, dec = self.coords

        if exposure is None:
            exposure = Exposure.createMockExposure(
                startTime=startTime, expTime=expTime, ra=ra, dec=dec,
                silent=silent, **kwargs)

        if exposure.isValid(silent=silent)[0] is False:
            if not silent:
                log.debug('mock exposure is invalid. Removing it.')
            return False

        validSet = logic.getValidSet(exposure, self)

        if validSet is not None:
            validSet.totoroExposures.append(exposure)
        else:
            self.sets.append(TotoroSet.Set.fromExposures([exposure]))

        return exposure

    def getIncompleteSet(self, startTime, expTime):
        """Returns incomplete sets that are valid for an exposure starting at
        JD=startTime."""

        LST0 = site.localSiderialTime(startTime)
        LST1 = (LST0 + expTime / 60.) % 24.

        incompleteSets = [set for set in self.sets if not set.complete]

        if len(incompleteSets) == 0:
            return None

        for set in incompleteSets:

            lstRange = set.getLSTRange()
            intersectionLength = logic.getIntervalIntersectionLength(
                (LST0, LST1), lstRange, wrapAt=24)

            if intersectionLength * 3600 > expTime:
                return set

        return None

    def getLastExposure(self):
        """Returns the last exposure taken."""

        exposures = self.getValidExposures()

        startTime = [exp.start_time for exp in exposures]
        order = np.argsort(startTime)

        return exposures[order[-1]]

    @property
    def priority(self):
        if not self.isMock:
            return self.plate_pointings[0].priority
        else:
            return 5

    @property
    def isPlugged(self):
        try:
            self.getActiveCartNumber()
            return True
        except:
            return False

    def getAPOcomplete(self):
        return utils.getAPOcomplete(self)

    def getMangaTileID(self):
        if len(self.design.inputs) == 0:
            return None
        else:
            for input in self.design.inputs:
                basename = os.path.basename(input.filepath)
                if 'mangaScience' not in basename:
                    continue
                noExt = os.path.splitext(basename)[0]
                return int(noExt.split('_')[1])
            return None

    def getAltitude(self, LST=None):
        """Returns the altitude of the plate at a certain LST."""

        if LST is None:
            LST = site.localSiderealTime()

        HA = (LST * 15. - self.ra) % 360.

        sinAlt = (np.sin(np.deg2rad(self.dec)) *
                  np.sin(np.deg2rad(site.latitude)) +
                  np.cos(np.deg2rad(self.dec)) *
                  np.cos(np.deg2rad(site.latitude)) *
                  np.cos(np.deg2rad(HA)))

        return np.rad2deg(np.arcsin(sinAlt))
