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
from Totoro.db import getConnectionFull
from Totoro import exceptions as TotoroExceptions
from Totoro import log, config, dustMap, site
from Totoro import utils
from Totoro.scheduler import observingPlan
from Totoro.dbclasses import Set as TotoroSet
from Totoro.dbclasses import Exposure as TotoroExposure
from Totoro.dbclasses import plate_utils as plateUtils
import warnings
from astropy import time
import numpy as np
from copy import deepcopy
from sqlalchemy import or_, and_
from sqlalchemy.orm.exc import NoResultFound


__all__ = ['getPlugged', 'getAtAPO', 'getAll', 'getComplete', 'Plate',
           'fromPlateID']


def getPlugged(**kwargs):

    kwargs.setdefault('fullCheck', False)
    kwargs.setdefault('updateSets', False)

    __, Session, plateDB, __ = getConnectionFull()
    session = Session()

    with session.begin():
        activePluggings = session.query(
            plateDB.ActivePlugging).join(
                plateDB.Plugging,
                plateDB.Plate,
                plateDB.PlateToSurvey,
                plateDB.Survey,
                plateDB.SurveyMode).filter(
                    plateDB.Survey.label == 'MaNGA',
                    plateDB.SurveyMode.label.ilike('MaNGA%')).order_by(
                        plateDB.Plate.plate_id).all()

    plates = [actPlug.plugging.plate for actPlug in activePluggings]

    return Plates(plates, **kwargs)


def getAtAPO(onlyIncomplete=False, onlyMarked=False,
             rejectLowPriority=False, raRange=None,
             rejectSpecial=True, **kwargs):
    """Gets plates at APO with various conditions."""

    kwargs.setdefault('fullCheck', False)
    kwargs.setdefault('updateSets', False)

    minimumPriority = config['plugger']['noPlugPriority']

    __, Session, plateDB, __ = getConnectionFull()
    session = Session()

    with session.begin():
        plates = session.query(plateDB.Plate).join(
            plateDB.PlateToSurvey, plateDB.Survey, plateDB.SurveyMode).filter(
                plateDB.Survey.label == 'MaNGA',
                plateDB.SurveyMode.label == 'MaNGA dither'
        ).join(plateDB.PlateLocation).filter(
            plateDB.PlateLocation.label == 'APO').join(
                plateDB.PlatePointing, plateDB.Pointing)

        _checkNumberPlatesAtAPO(plates)

        if onlyMarked is False:
            plates = plates
        else:
            plates = plates.join(plateDB.PlateToPlateStatus,
                                 plateDB.PlateStatus).filter(
                plateDB.PlateStatus.label == 'Accepted')

        if raRange is not None:
            if raRange.size == 2:
                plates = plates.filter(
                    plateDB.Pointing.center_ra >= raRange[0],
                    plateDB.Pointing.center_ra <= raRange[1])
            elif raRange.size == 4:
                plates = plates.filter(
                    or_(and_(plateDB.Pointing.center_ra >= raRange[0][0],
                             plateDB.Pointing.center_ra <= raRange[0][1]),
                        and_(plateDB.Pointing.center_ra >= raRange[1][0],
                             plateDB.Pointing.center_ra <= raRange[1][1])))
            else:
                warnings.warn('unrecognised format for raRange',
                              TotoroExceptions.TotoroUserWarning)

        if rejectLowPriority:
            plates = plates.filter(
                plateDB.PlatePointing.priority > minimumPriority)

        plates = plates.order_by(plateDB.Plate.plate_id).all()

        validPlates = []
        if rejectSpecial:
            for plate in plates:
                if plate.mangadbPlate is None:
                    continue
                if (plate.mangadbPlate.all_sky_plate is True or
                        plate.mangadbPlate.commissioning_plate is True or
                        plate.mangadbPlate.neverobserve is True):
                    continue
                validPlates.append(plate)
        else:
            validPlates = plates

    if onlyIncomplete:
        return _getIncomplete(validPlates, **kwargs)
    else:
        return Plates(validPlates, **kwargs)


def getAll(onlyIncomplete=False, **kwargs):

    kwargs.setdefault('fullCheck', False)
    kwargs.setdefault('updateSets', False)

    __, Session, plateDB, __ = getConnectionFull()
    session = Session()

    with session.begin():
        plates = session.query(plateDB.Plate).join(
            plateDB.PlateToSurvey, plateDB.Survey, plateDB.SurveyMode
        ).filter(plateDB.Survey.label == 'MaNGA',
                 plateDB.SurveyMode.label == 'MaNGA dither').order_by(
                     plateDB.Plate.plate_id)

    platesAtAPO = plates.join(plateDB.PlateLocation).filter(
        plateDB.PlateLocation.label == 'APO')
    _checkNumberPlatesAtAPO(platesAtAPO)

    plates = plates.order_by(plateDB.Plate.plate_id).all()

    if onlyIncomplete:
        return _getIncomplete(plates, **kwargs)
    else:
        return Plates(plates, **kwargs)


def _getIncomplete(plates, **kwargs):

    totoroPlates = Plates(plates, **kwargs)

    incompletePlates = []
    for totoroPlate in totoroPlates:
        if not totoroPlate.isComplete:
            incompletePlates.append(totoroPlate)

    return incompletePlates


def getComplete(**kwargs):
    """Retrieves only complete files by looking at the plugging status
    (faster than using getAll and then selecting the complete plates)."""

    kwargs.setdefault('fullCheck', False)
    kwargs.setdefault('updateSets', False)

    db, Session, plateDB, mangaDB = getConnectionFull()
    session = Session()

    with session.begin():
        plates = session.query(plateDB.Plate).join(
            plateDB.PlateToSurvey, plateDB.Survey, plateDB.SurveyMode,
            plateDB.Plugging, plateDB.PluggingStatus).filter(
                plateDB.Survey.label == 'MaNGA',
                plateDB.SurveyMode.label == 'MaNGA dither',
                plateDB.PluggingStatus.label.in_(['Good', 'Overridden Good'])
        ).order_by(plateDB.Plate.plate_id).all()

    return Plates(plates, **kwargs)


def _checkNumberPlatesAtAPO(plates):
    """Raises a warning if the number of plates at APO is reaching the limit
    we can store there."""

    nPlatesAPO = config['numberPlatesAllowedAtAPO']
    nPlates = plates.count()
    if nPlates >= 0.9 * nPlatesAPO:
        warnings.warn('MaNGA has {0} plates at APO when the maximum is {1}'
                      .format(nPlates, nPlatesAPO),
                      TotoroExceptions.TotoroUserWarning)

    return


class Plates(list):

    def __init__(self, inp, format='pk', **kwargs):

        __, __, plateDB, __ = getConnectionFull()

        if all([isinstance(ii, Plate) for ii in inp]):
            list.__init__(self, inp)
        elif all([isinstance(ii, plateDB.Plate) for ii in inp]):
            list.__init__(self, [Plate(ii, **kwargs) for ii in inp])
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
        return getAll(**kwargs)

    @staticmethod
    def getAtAPO(**kwargs):
        """For backards compatibility."""
        return getAtAPO(**kwargs)


def fromPlateID(plateid, **kwargs):
    """Convenience function that returns a `Plate` instance from a plate_id."""
    return Plate(input=plateid, format='plate_id', **kwargs)


class Plate(object):

    _instances = {}

    def __new__(cls, input=None, format='pk', **kwargs):

        # Checks if an open transaction already exists.
        utils.checkOpenSession()

        me = object.__new__(cls)

        me.db, Session, plateDB, __ = getConnectionFull()
        me.session = Session()

        # Initialises DB object
        if input is None:
            me._dbObject = plateDB.Plate()
            return me
        else:
            if isinstance(input, plateDB.Plate):
                me._dbObject = input
            else:
                me._dbObject = me._initFromData(input, format=format)

        # If the DB object already exists in the library of Totoro.Plate
        # instances, returns it. Otherwise, records it and returns the new
        # object.
        if me._dbObject in cls._instances:
            return cls._instances[me._dbObject]
        else:
            cls._instances[me._dbObject] = me
            return me

    def __init__(self, input=None, format='pk', mock=False,
                 updateSets=True, mjd=None, fullCheck=True,
                 manga_tileid=None, **kwargs):
        """A custom class based on plateDB.Plate."""

        self.__dbAttributes__ = self._dbObject.__mapper__.attrs

        self._complete = None
        self._drilled = None
        self._priority = None
        self.isMock = mock
        self._kwargs = kwargs
        self.mjd = mjd
        self._manga_tileid = manga_tileid

        # Date (JD) at which the plate will arrive at APO. Used for field
        # selection. If None, it assumes the plate is already at APO.
        self._dateAtAPO = None

        if 'dust' in kwargs:
            self.dust = kwargs['dust']
        else:
            self.dust = None if dustMap is None else dustMap(self.ra, self.dec)

        self.mlhalimit = utils.mlhalimit(self.dec)

        if not self.isMock:
            self.sets = [TotoroSet(set, **kwargs)
                         for set in self.getMangaDBSets()]

            self.checkPlate(full=fullCheck)

            if updateSets:
                self.updatePlate(**kwargs)

        else:
            self.sets = []

    def __repr__(self):
        return ('<Totoro Plate (plate_id={0}, manga_tileid={1}, '
                'completion={2:.2f})>'
                .format(self.plate_id, self.manga_tileid,
                        self.getPlateCompletion()))

    def _initFromData(self, input, format):
        """Init a new Plate instance from a DB query."""

        with self.session.begin():
            try:
                plate = self.session.query(self.db.plateDB.Plate).filter(
                    eval('self.db.plateDB.Plate.{0} == {1}'
                         .format(format, input))).one()
            except NoResultFound:
                raise TotoroExceptions.TotoroError('no plate found for input '
                                                   '{0}={1}'
                                                   .format(format, input))

        return plate

    def __getattr__(self, name):
        """Custom getattr method that first looks into the DB object."""

        if name in object.__getattribute__(self, '__dbAttributes__'):
            return getattr(self._dbObject, name)
        else:
            return object.__getattribute__(self, name)

    def __setattr__(self, name, value):
        """Custom setattr method that first looks into the DB object."""

        if hasattr(self, '_dbObject') and hasattr(self._dbObject, name):
            setattr(self._dbObject, name, value)
        else:
            super(Plate, self).__setattr__(name, value)

    @classmethod
    def fromSets(cls, sets, **kwargs):

        newPlate = cls(None, mock=True, **kwargs)
        newPlate.sets = sets

        log.debug('created mock plate from sets pk={0}'.format(
            ', '.join(map(str, sets))))

        return newPlate

    @classmethod
    def createMockPlate(cls, ra=None, dec=None, **kwargs):

        if ra is None or dec is None:
            raise TotoroExceptions.TotoroError('ra and dec must be specified')

        mockPlate = cls(mock=True, ra=ra, dec=dec, **kwargs)
        mockPlate.isMock = True

        log.debug('created mock plate with ra={0:.3f} and dec={0:.3f}'.format(
                  mockPlate.coords[0], mockPlate.coords[1]))

        return mockPlate

    def addExposure(self, exposure):
        return plateUtils.addExposure(exposure, self)

    def getMangaDBSets(self):
        """Returns a list of mangaDB.ModelClasses.Set instances with all the
        sets matching exposures in this plate. Removes duplicates."""

        plateDB = self.db.plateDB
        mangaDB = self.db.mangaDB

        # Finds the sets for this plate
        sets = self.session.query(mangaDB.Set).join(
            mangaDB.Exposure, plateDB.Exposure, plateDB.Observation,
            plateDB.PlatePointing, plateDB.Plate).filter(
                plateDB.Plate.pk == self.pk).all()

        sets = sorted(sets, key=lambda set: set.pk)

        return sets

    def checkPlate(self, full=True):
        """Does some sanity checks for the current plate. If full=True,
        it checks that all the science exposures have a mangaDB counterpart.
        This can be disabled to save time during plate bulk load."""

        # Checks that we are not inside an open session
        utils.checkOpenSession()

        if not hasattr(self, 'pk') or not hasattr(self, 'plate_id'):
            raise AttributeError('Plate instance has no pk or plate_id.')

        if not self.isMaNGA:
            raise TotoroExceptions.NoMangaPlate('this is not a MaNGA plate!')

        if full:
            nMaNGAExposures = len(self.getMangadbExposures())
            nScienceExposures = len(self.getScienceExposures())
            if nMaNGAExposures != nScienceExposures:
                warnings.warn('plate_id={1}: {0} plateDB.Exposures found '
                              'but only {2} mangaDB.Exposures'.format(
                                  nScienceExposures, self.plate_id,
                                  nMaNGAExposures),
                              TotoroExceptions.NoMangaExposure)

    @property
    def isMaNGA(self):

        for survey in self.surveys:
            if survey.label == 'MaNGA':
                return True

        return False

    def updatePlate(self, force=False, **kwargs):

        # Checks if some sets have incompletely reduced exposures. If so,
        # removes their set_pk so that they are completely procesed again.
        # with session.begin():
        #     for ss in self.sets:
        #         for exp in ss.totoroExposures:
        #             if exp._mangaExposure.status is None:
        #                 if not exp.isMock:
        #                     exp.mangadbExposure[0].set_pk = None
        #                 ss.totoroExposures.remove(exp)
        #                 ss.set_status_pk = None

        #         # Checks if the sets is now empty, if so, removes it
        #         # from the object and from the DB
        #         if len(ss.totoroExposures) == 0:
        #             if not ss.isMock:
        #                 setDB = session.query(mangaDB.Set).get(ss.pk)
        #                 session.delete(setDB)
        #             self.sets.remove(ss)

        # if self.isComplete:
        #     if not force:
        #         log.debug('plate_id={0} is marked complete. '
        #                   'Not updating sets.'.format(self.plate_id))
        #         return False
        #     else:
        #         log.info('plate_id={0} is marked complete but force=True.'
        #                  .format(self.plate_id))

        if (self.currentSurveyMode is None or
                self.currentSurveyMode.label != 'MaNGA dither'):
            log.debug('plate_id={0} has surveyMode which is not MaNGA dither. '
                      'Not updating sets.')
            return False

        result = plateUtils.updatePlate(self, **kwargs)

        return result

    def rearrangeSets(self, LST=None, mode='optimal', scope='all', **kwargs):
        result = plateUtils.rearrangeSets(self, LST=LST, scope=scope,
                                          mode=mode, **kwargs)

        return result

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
                              TotoroExceptions.MultiplePlatePointings)

            return np.array(
                [self.plate_pointings[0].pointing.center_ra,
                 self.plate_pointings[0].pointing.center_dec], np.float)

    @property
    def isComplete(self):
        if self._complete is not None:
            return self._complete
        else:
            status = utils.isPlateComplete(self)
            # If the plate is mock, caches the status.
            if status is True and self.isMock:
                self._complete = status
            return status

    def copy(self):
        return deepcopy(self)

    def getPlateCompletion(self, includeIncompleteSets=False, useMock=True):

        # totalSN = self.getCumulatedSN2(
        #     includeIncomplete=includeIncompleteSets, useMock=useMock)

        validStatuses = ['Good', 'Excellent', 'Override Good']
        if includeIncompleteSets:
            validStatuses.append('Incomplete')

        sets = []
        for ss in self.sets:
            if not useMock and ss.isMock:
                continue
            if ss.getStatus()[0] in validStatuses:
                sets.append(ss)

        exposureSN2 = np.array([exp.getSN2Array() for ss in sets
                                for exp in ss.totoroExposures])

        if exposureSN2.shape[0] == 0:
            return 0

        blueSN2 = np.nansum(np.nanmean(exposureSN2[:, [0, 1]], axis=1))
        redSN2 = np.nansum(np.nanmean(exposureSN2[:, [2, 3]], axis=1))

        if np.isnan(blueSN2) or np.isnan(redSN2):
            return 0
        else:
            return np.min([blueSN2 / config['SN2thresholds']['plateBlue'],
                          redSN2 / config['SN2thresholds']['plateRed']])

    def getSN2Array(self, **kwargs):
        """Same as getCumulatedSN2. Added for consistency with Set and
        Exposure."""
        return self.getCumulatedSN2(**kwargs)

    def getCumulatedSN2(self, includeIncomplete=False, useMock=True,
                        **kwargs):
        """Returns the cumulated SN2 from all valid sets.
        If includeIncomplete=True, incomplete sets SN2 are also included.
        If useMock=False, mock sets are ignored.
        """

        validStatuses = ['Good', 'Excellent', 'Override Good']
        if includeIncomplete:
            validStatuses.append('Incomplete')

        validSets = []
        for set in self.sets:
            # Creates a mock set with the appropriate exposures.
            if useMock:
                mockSet = set
            else:
                if set.isMock:
                    continue
                mockSet = TotoroSet.fromExposures(
                    [exp for exp in set.totoroExposures if not exp.isMock])
            if mockSet.getQuality()[0] in validStatuses:
                validSets.append(mockSet)

        if len(validSets) == 0:
            return np.array([0.0, 0.0, 0.0, 0.0])
        else:
            return np.nansum([set.getSN2Array(useMock=useMock)
                              for set in validSets], axis=0)

    def getActiveCartNumber(self):
        """Returns the cart number of the active plugging. Raises an error if
        no active pluggings are found."""

        for plugging in self.pluggings:
            if len(plugging.activePlugging) > 0:
                return int(plugging.cartridge.number)

        raise TotoroExceptions.PlateNotPlugged(
            'plate_id={0} is not currently plugged'.format(self.plate_id))

    def getActivePlugging(self):
        """Returns the active plugging or None if none found."""

        if self.isMock:
            return None

        for plugging in self.pluggings:
            if len(plugging.activePlugging) > 0:
                return plugging

        return None

    def getMangadbExposures(self):
        """Returns a list of mangaDB.Exposure objects with all the exposures
        for this plate."""

        scienceExps = self.getScienceExposures()
        mangaExposures = []
        for exp in scienceExps:
            if exp.mangadbExposure is None or len(exp.mangadbExposure) == 0:
                raise TotoroExceptions.TotoroError(
                    'platedb exposure_no={0} has no mangadb counterpart'
                    .format(exp.exposure_no))
            else:
                mangaExposures.append(exp.mangadbExposure[0])

        return mangaExposures

    def getScienceExposures(self):
        """Returns a list of all plateDB.Exposure science exposures for
        the plate."""

        scienceExps = []

        plateDB = self.db.plateDB

        for plugging in self.pluggings:
            exposures = self.session.query(plateDB.Exposure).join(
                plateDB.Observation, plateDB.ExposureFlavor).filter(
                    plateDB.Observation.plugging_pk == plugging.pk).filter(
                        plateDB.ExposureFlavor.label == 'Science').all()
            scienceExps += exposures

        return scienceExps

    def getTotoroExposures(self, onlySets=False):
        """Returns a list of Totoro dbclasses.Exposure instances for this
        Plate insance."""

        totoroExposures = []
        for ss in self.sets:
            totoroExposures += ss.totoroExposures

        if onlySets or self.isMock:
            return totoroExposures

        # Some exposures not in sets may be missing.
        allExposures = self.getMangadbExposures()
        totoroExposures += [TotoroExposure(exp.pk, parent='mangadb')
                            for exp in allExposures if exp.set_pk is None]

        return totoroExposures

    def getValidSets(self, includeIncomplete=False):

        validSets = []
        for set in self.sets:
            quality = set.getQuality()[0]
            if quality in ['Good', 'Excellent', 'Override Good']:
                validSets.append(set)
            elif includeIncomplete and quality == 'Incomplete':
                validSets.append(set)

        return validSets

    def getValidExposures(self, **kwargs):
        """Returns all valid exposures, even if they belong to an incomplete
        or bad set."""

        validExposures = []

        for set in self.sets:
            for exp in set.totoroExposures:
                if exp.isValid(**kwargs)[0] is True:
                    validExposures.append(exp)

        return validExposures

    def getHARange(self, intersect=False, mjd=None, **kwargs):

        ha0 = -self.mlhalimit
        ha1 = self.mlhalimit

        haRange = np.array([ha0, ha1]) % 360.
        haRange[haRange > 180.] -= 360.

        if intersect is False:
            return haRange

        if mjd is None:
            mjd = int(np.round(time.Time.now().mjd))
            log.debug('Plate.getHARange: MJD set to {0:d}'.format(mjd))
        else:
            mjd = int(mjd)

        jdRange = observingPlan.getMJD(mjd)

        if jdRange is None:
            warnings.warn('no observing block found for MJD={0:d}. '
                          'Observing windows will not be contrained.'
                          .format(mjd), TotoroExceptions.NoObservingBlock)
            return haRange

        observingRangeLST = np.array(map(site.localSiderealTime, jdRange))
        observingRangeHA = (observingRangeLST * 15 - self.ra) % 360.

        haRange = utils.getIntervalIntersection(haRange, observingRangeHA)
        haRange[haRange > 180] -= 360.

        return haRange

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
                        expTime=None, silent=False, rearrange=True, **kwargs):
        """Creates a mock expusure in the best possible way."""

        ra, dec = self.coords

        if exposure is None:
            exposure = TotoroExposure.createMockExposure(
                startTime=startTime, expTime=expTime, ra=ra, dec=dec,
                dust=self.dust, **kwargs)

        validSet = plateUtils.getOptimalSet(self, exposure)

        if validSet is not None:
            # Gets a valid dither position for the exposure in this set.
            dither = plateUtils.getValidDither(validSet)
            exposure.ditherPosition = dither
            newSet = False

        else:
            # This is a new set. We give the exposure a random dither position.
            exposure.ditherPosition = 'N'
            validSet = TotoroSet(mock=True, **kwargs)
            newSet = True

        if exposure.isValid()[0] is False:
            if not silent:
                log.debug('mock exposure is invalid. Removing it.')
            return False

        if newSet:
            self.sets.append(validSet)
        validSet.totoroExposures.append(exposure)

        if rearrange:
            # We only run this rearrangement if the list of exposures in
            # incomplete sets is small
            nExpIncompleteSets = len(self.getExposuresInIncompleteSets())
            if nExpIncompleteSets <= 5:
                self.rearrangeSets(mode='optimal', scope='incomplete',
                                   silent=True)
            else:
                warnings.warn('plate={0}: skipping incomplete set '
                              'rearrangement because plate has > 5 exposures '
                              'in incomplete sets',
                              TotoroExceptions.TotoroUserWarning)

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
            intersectionLength = utils.getIntervalIntersectionLength(
                (LST0, LST1), lstRange, wrapAt=24)

            if intersectionLength * 3600 > expTime:
                return set

        return None

    def getExposuresInIncompleteSets(self):
        """Returns a list of exposures in incomplete sets."""

        incompleteSetExposures = []

        incompleteSets = [set for set in self.sets if not set.complete]
        for ss in incompleteSets:
            incompleteSetExposures += ss.totoroExposures

        return incompleteSetExposures

    def getLastExposure(self):
        """Returns the last exposure taken."""

        exposures = self.getValidExposures()

        startTime = [exp.start_time for exp in exposures]
        order = np.argsort(startTime)

        return exposures[order[-1]]

    def getLastSet(self):
        """Returns the set containing the last exposure."""

        lastExposure = self.getLastExposure()

        for ss in self.sets:
            if lastExposure in ss.totoroExposures:
                return ss

    @property
    def priority(self):
        if not self.isMock and self._priority is None:
            return self.plate_pointings[0].priority
        else:
            if self._priority is None:
                return config['defaultPriority']
            else:
                return self._priority

    @priority.setter
    def priority(self, value):
        if not self.isMock:
            with self.session.begin():
                self.plate_pointings[0].priority = value
        else:
            self._priority = value

    @property
    def isPlugged(self):
        try:
            self.getActiveCartNumber()
            return True
        except:
            return False

    def getAPOcomplete(self, **kwargs):
        """Retuns the APOcomplete dictionary for this plate."""
        return utils.getAPOcomplete(self, **kwargs)

    def createAPOcompleteFile(self, path=None):
        """Creates an apocomp file for the current plate."""

        return utils.createAPOcompleteFile(self.getAPOcomplete(), path=path)

    @property
    def manga_tileid(self):
        """Returns manga_tileid for this plate."""
        if self._manga_tileid is not None:
            return self._manga_tileid

        if self.mangadbPlate is None:
            return None
        else:
            self._manga_tileid = self.mangadbPlate.manga_tileid
            return self._manga_tileid

    @manga_tileid.setter
    def manga_tileid(self, value):
        """Sets manga_tileid for this plate. Not saved to the DB."""
        self._manga_tileid = value

    def getMangaTileID(self):
        """Returns manga_tiled. Backwards compatibility."""
        return self.manga_tileid

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

    def getLSTRangeAboveAltitude(self, altitude):
        """Returns the LST range in which the plate is above ``altitude`` [deg]
        """

        altRad = np.deg2rad(altitude)
        latRad = np.deg2rad(site.latitude)
        decRad = np.deg2rad(self.dec)

        xx = ((np.sin(altRad) - np.sin(latRad) * np.sin(decRad)) /
              (np.cos(latRad) * np.cos(decRad)))

        if xx > 1 or xx < -1:
            return np.array([0, 0])

        HA = np.rad2deg(np.arccos(xx))

        LSTRange = [(-HA + self.ra) % 360. / 15., (HA + self.ra) % 360. / 15.]

        return np.array(LSTRange)

    def getNumberPermutations(self):
        """Returns the number of permutations for a brute-force set
        rearrangement."""

        validExposures = self.getValidExposures()
        ditherPositions = [exp.ditherPosition for exp in validExposures]
        return plateUtils.getNumberPermutations(ditherPositions)

    def getPermutations(self):
        """Returns the permutations for a brute-force set rearrangement."""

        validExposures = self.getValidExposures()
        ditherPositions = [exp.ditherPosition for exp in validExposures]
        return plateUtils.calculatePermutations(ditherPositions)

    def getLocation(self):
        """Returns the location label of the plate."""

        if self.location is None:
            return None
        return self.location.label

    def getOrphaned(self, includeUnplugged=True, useMock=True):
        """Returns orphaned exposures in sets. Exposures must be valid
        and belong to a valid set. If `useMock=False`, only real sets will
        be considered. If `includeUnplugged=True`, orphaned exposures in
        unplugged sets will be included."""

        orphaned = []
        validStatuses = ['Incomplete']
        if includeUnplugged:
            validStatuses.append('Unplugged')

        for set in self.sets:
            # Creates a mock set with the appropriate exposures
            if not useMock:
                if set.isMock:
                    continue
                mockSet = TotoroSet.fromExposures(
                    [exp for exp in set.totoroExposures if not exp.isMock])
            else:
                mockSet = set

            # Checks if set is incomplete. If so, the exposures are orphaned.
            if mockSet.getQuality()[0] in validStatuses:
                orphaned += mockSet.totoroExposures

        return orphaned

    def getMockExposures(self):
        """Returns mock exposues in this plate."""

        result = []
        for ss in self.sets:
            for exp in ss.totoroExposures:
                if exp.isMock:
                    result.append(exp)

        return result

    def hasIncompleteSets(self):
        """Returns True if the plate has at least one incomplete set."""

        for ss in self.sets:
            if ss.getStatus()[0] == 'Incomplete':
                return True

        return False

    @property
    def drilled(self):
        """Property to record if a plate has already been drilled."""
        if self._drilled is not None:
            return self._drilled
        else:
            if self.plate_id is not None:
                return True
            else:
                return False

    @drilled.setter
    def drilled(self, value):
        self._drilled = value

    @property
    def dateAtAPO(self):
        return self._dateAtAPO

    @dateAtAPO.setter
    def dateAtAPO(self, value):
        assert value >= 0, 'value must be a JD.'
        self._dateAtAPO = value
