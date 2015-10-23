#!/usr/bin/env python
# encoding: utf-8
"""
set.py

Created by José Sánchez-Gallego on 18 Mar 2014.
Licensed under a 3-clause BSD license.

Revision history:
    18 Mar 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from exposure import Exposure
from Totoro import TotoroDBConnection
from Totoro import log, config, site
from Totoro import exceptions
from Totoro import utils
# from Totoro import scheduler
import numpy as np
import warnings
from copy import copy
from astropy import time


db = TotoroDBConnection()
plateDB = db.plateDB
mangaDB = db.mangaDB
session = db.Session()


__all__ = ['Set']


def getPlateSets(inp, format='plate_id', **kwargs):

    with session.begin():
        sets = session.query(mangaDB.Set).join(
            mangaDB.Exposure,
            plateDB.Exposure,
            plateDB.Observation,
            plateDB.PlatePointing,
            plateDB.Plate).filter(
                eval('plateDB.Plate.{0} == {1}'.format(format, inp))).all()

    return [Set(set, **kwargs) for set in sets]


class Set(mangaDB.Set):

    def __new__(cls, input=None, format='pk', *args, **kwargs):

        if input is None:
            ss = mangaDB.Set.__new__(cls)
            super(Set, ss).__init__(**kwargs)
            return ss

        base = cls.__bases__[0]

        if isinstance(input, base):
            instance = input
        else:
            with session.begin():
                instance = session.query(base).filter(
                    eval('mangaDB.Set.{0} == {1}'.format(format, input))).one()

        instance.__class__ = cls

        return instance

    def __init__(self, inp=None, format='pk', mock=False, mjd=None,
                 *args, **kwargs):

        self._status = None

        self.isMock = mock
        if inp is None:
            self.isMock = True

        self._kwargs = kwargs
        self.mjd = mjd

        if not self.isMock:
            self.totoroExposures = self.loadExposures()
        else:
            self.totoroExposures = []

        # Checks that the set has exposures.
        if not self.isMock:
            self._checkHasExposures()

    def __repr__(self):
        return '<Totoro Set (pk={0}, status={1})>'.format(
            self.pk, self.getQuality(flag=False)[0])

    def update(self, **kwargs):
        """Reloads the set."""

        newSelf = Set(self.pk, fromat='pk', mock=self.isMock,
                      mjd=self.mjd, **self._kwargs)
        self = newSelf

        log.debug('Set pk={0} has been reloaded'.format(
                  self.pk))

    @staticmethod
    def fromExposures(exposures, **kwargs):
        """Creates a mock set for a list of Totoro.Exposures."""

        newSet = Set(mock=True, **kwargs)

        if isinstance(exposures, Exposure):
            newSet.totoroExposures = [exposures]
        else:
            for exp in exposures:
                assert isinstance(exp, Exposure), \
                    '{0} not an instance of Totoro.Exposure'.format(exp)
            newSet.totoroExposures = list(exposures)

        return newSet

    def loadExposures(self):

        return [Exposure(mangaExp) for mangaExp in self.exposures]

    def _checkHasExposures(self):
        if len(self.totoroExposures) == 0:
            raise exceptions.EmptySet(
                'set pk={0} has no exposures.'.format(self.pk))

    @classmethod
    def createMockSet(cls, ra=None, dec=None, **kwargs):

        if ra is None or dec is None:
            raise exceptions.TotoroError('ra and dec must be specified')

        newSet = Set(mock=True, ra=ra, dec=dec, **kwargs)

        return newSet

    def addMockExposure(self, **kwargs):

        if self.complete is True:
            raise exceptions.TotoroError(
                'set is complete; not exposures can be added.')

        if 'ditherPosition' not in kwargs or kwargs['ditherPosition'] is None:
            kwargs['ditherPosition'] = self.getMissingDitherPositions()[0]

        newExposure = Exposure.createMockExposure(**kwargs)
        self.totoroExposures.append(newExposure)

    @property
    def plate(self):
        """Returns the plateDB.Plate object for this set."""

        self._checkHasExposures()
        return (self.exposures[0].platedbExposure.observation.
                plate_pointing.plate)

    @property
    def ra(self):
        return self.getCoordinates()[0]

    @property
    def dec(self):
        return self.getCoordinates()[1]

    def getCoordinates(self):

        if 'ra' in self._kwargs and 'dec' in self._kwargs:
            if (self._kwargs['ra'] is not None and
                    self._kwargs['dec'] is not None):
                return np.array(
                    [self._kwargs['ra'], self._kwargs['dec']], np.float)
        else:
            self._checkHasExposures()
            return self.totoroExposures[0].getCoordinates()

    def getHA(self, **kwargs):
        """Returns the HA interval of the exposures in the set. If midPoint is
        set, the middle point of the exposures is used for the calculation."""

        exposures = self.getValidExposures(**kwargs)

        if len(exposures) == 0:
            plateHALimit = utils.mlhalimit(self.dec)
            HA = np.array([-plateHALimit, plateHALimit])
        elif len(exposures) >= 1:
            expHAs = np.array([exp.getHA() for exp in exposures])
            HA = np.array(utils.getMinMaxIntervalSequence(expHAs))

        HA[HA > 180] -= 360
        return HA

    def getHARange(self, intersect=False, mjd=None,
                   maxHARange=config['set']['maxHARange'], **kwargs):
        """Returns the HA limits to add more exposures to the set."""

        ha = self.getHA()
        haRange = np.array([np.max(ha) - maxHARange,
                            np.min(ha) + maxHARange]) % 360.

        plateHALimit = utils.mlhalimit(self.dec)

        haRangePlate = utils.getIntervalIntersection(
            haRange, np.array([-plateHALimit, plateHALimit]), wrapAt=360.)
        haRangePlate[haRangePlate > 180] -= 360

        # if intersect is False:
        return haRangePlate

        # if mjd is None:
        #     mjd = int(np.round(time.Time.now().mjd))
        #     log.debug('Set.getHARange: MJD set to {0:d}'.format(mjd))
        # else:
        #     mjd = int(mjd)

        # jdRange = scheduler.observingPlan.getMJD(mjd)
        # mangaLSTRange = np.array(site.localSiderealTime(jdRange))

        # mangaHARange = (mangaLSTRange * 15. - self.ra) % 360.

        # haRange = utils.getIntervalIntersection(haRangePlate, mangaHARange,
        #                                         wrapAt=360.)
        # if haRange is False:
        #     return False
        # else:
        #     return haRange

    def getDitherPositions(self):
        """Returns a list of dither positions in the set."""

        return [exp.ditherPosition for exp in self.totoroExposures]

    def getMissingDitherPositions(self):
        """Returns a list of missing dither positions."""

        setDitherPositions = copy(config['set']['ditherPositions'])
        for dPos in self.getDitherPositions():
            if dPos.upper() in setDitherPositions:
                setDitherPositions.remove(dPos.upper())

        return setDitherPositions

    def getSN2Array(self, **kwargs):
        """Returns an array with the cumulated SN2 of the valid exposures in
        the set. The return format is [b1SN2, b2SN2, r1SN2, r2SN2]."""

        if len(self.totoroExposures) == 0:
            return np.array([0.0, 0.0, 0.0, 0.0])
        else:
            return np.nansum([exp.getSN2Array()
                              for exp in self.totoroExposures], axis=0)

    def getSN2Range(self):
        """Returns the SN2 range in which new exposures may be taken."""

        maxSN2Factor = config['set']['maxSN2Factor']

        sn2 = np.array([exp.getSN2Array() for exp in self.totoroExposures])
        sn2Average = np.array(
            [(np.nanmean(ss[0:2]), np.nanmean(ss[2:4])) for ss in sn2])

        minSN2Blue = np.max(sn2Average[:, 0]) / maxSN2Factor
        maxSN2Blue = np.min(sn2Average[:, 0]) * maxSN2Factor
        minSN2Red = np.max(sn2Average[:, 1]) / maxSN2Factor
        maxSN2Red = np.min(sn2Average[:, 1]) * maxSN2Factor

        if minSN2Red < config['SN2thresholds']['exposureRed']:
            minSN2Red = config['SN2thresholds']['exposureRed']
        if minSN2Blue < config['SN2thresholds']['exposureBlue']:
            minSN2Blue = config['SN2thresholds']['exposureBlue']

        return np.array([[minSN2Blue, maxSN2Blue], [minSN2Red, maxSN2Red]])

    def getSeeingRange(self):
        """Returns the seeing range in which new exposures may be taken."""

        maxSeeingRange = config['set']['maxSeeingRange']
        goodSeeing = config['set']['goodSeeing']
        maxSeeing = config['exposure']['maxSeeing']

        nExposuresMissing = 3 - len(self.totoroExposures)
        seeings = np.array([exp.seeing for exp in self.totoroExposures])

        maxSeeingGood = (3 * goodSeeing - np.sum(seeings)) / nExposuresMissing
        if maxSeeingGood < maxSeeing:
            maxSeeing = maxSeeingGood

        seeingRangeMin = np.max(seeings) - maxSeeingRange
        seeingRangeMax = np.min(seeings) + maxSeeingRange
        if seeingRangeMax > maxSeeing:
            seeingRangeMax = maxSeeing

        return np.array([seeingRangeMin, seeingRangeMax])

    def getQuality(self, **kwargs):
        """Alias for getStatus."""

        return self.getStatus(**kwargs)

    def getStatus(self, **kwargs):
        """Returns the status of the set."""

        if self._status is not None:
            return self._status

        status = checkSet(self, **kwargs)

        if self.isMock and status[0] != 'Incomplete':
            self._quality = status

        return status

    def getValidExposures(self, **kwargs):

        validExposures = []
        for exp in self.totoroExposures:
            if exp.isValid(**kwargs)[0] is True:
                validExposures.append(exp)

        return validExposures

    def getAverageSeeing(self, **kwargs):

        seeings = []
        for exp in self.totoroExposures:
            if exp.isValid(**kwargs)[0]:
                seeings.append(exp.seeing)

        return np.mean(seeings)

    def getLST(self):

        ha0, ha1 = self.getHA()

        lst0 = (ha0 + self.ra) % 360. / 15
        lst1 = (ha1 + self.ra) % 360. / 15

        return np.array([lst0, lst1])

    def getLSTRange(self, **kwargs):

        haRange = self.getHARange(**kwargs)

        if haRange is False:
            return False

        ha0, ha1 = haRange

        lst0 = (ha0 + self.ra) % 360. / 15
        lst1 = (ha1 + self.ra) % 360. / 15

        return np.array([lst0, lst1])

    def getUTVisibilityWindow(self, **kwargs):
        return self.getUTRange(**kwargs)

    def getUTRange(self, mjd=None, returnType='str', **kwargs):

        lstRange = self.getLSTRange(mjd=mjd, **kwargs)

        if lstRange is False:
            return False

        lst0, lst1 = lstRange

        if mjd is None:
            mjd = int(np.round(time.Time.now().mjd))

        date = time.Time(mjd, format='mjd', scale='tai')

        date0 = site.localSiderealTimeToDate(lst0, date=date)
        date1 = date0 + time.TimeDelta(
            (lst1 - lst0) % 24. * 3600., format='sec')

        if returnType == 'str':
            return ('{0:%H:%M}'.format(date0.datetime),
                    '{0:%H:%M}'.format(date1.datetime))
        else:
            return (date0.datetime, date1.datetime)

    @property
    def complete(self):
        if self.getQuality()[0] not in ['Incomplete', 'Bad']:
            return True
        else:
            return False


def flagSet(set, statusLabel, errorCode, flag=True, message=None,
            silent=False, **kwargs):
    """Helper function to log and flag sets."""

    if message is not None and not silent:
        log.debug(message)

    if flag:
        # Avoids calling setSetStatus if there is nothing to flag
        if (statusLabel in ['Incomplete', 'Unplugged', 'Bad'] and
                set.set_status_pk is None):
            pass
        else:
            setSetStatus(set, statusLabel)

    return (statusLabel, errorCode)


def setSetStatus(set, status):
    """ Sets the status of an exposure.

    Parameters
    ----------
    set : `Totoro.dbclasses.Set` or `mangaDB.Set` instance, or integer.
        The set that will receive the new status. Either a
        Totoro.dbclasses.Set, mangaDB.Set, or the pk of the mangaDB.Set.
    status : string
        The status to be set. It must be one of the values in
        mangaDB.SetStatus.label.

    Returns
    -------
    result : bool
        Returns True if the status has been set correctly.

    """

    # from Totoro.dbclasses import Set

    if isinstance(set, (Set, db.mangaDB.Set)):
        pk = set.pk
    else:
        pk = set

    with session.begin():
        try:
            queryStatus = session.query(db.mangaDB.SetStatus).filter(
                db.mangaDB.SetStatus.label == status).one()
            statusPK = queryStatus.pk
        except:
            # If the status is not found, we remove the status.
            statusPK = None

        ss = session.query(db.mangaDB.Set).get(pk)

        if ss.set_status_pk is not None and statusPK is None:
            warnings.warn('changing set pk={0} from status {1} to None'
                          .format(pk, ss.status.label),
                          exceptions.TotoroUserWarning)

        ss.set_status_pk = statusPK

    log.debug('mangaDB.Set.pk={0} set to {1}'.format(pk, status))

    return True


setErrorCodes = {
    0: 'no error.',
    1: 'one or more exposures are invalid.',
    2: 'HA range is greater than maximum allowed.',
    3: 'seeing values out of range',
    4: 'SN2 values out of range',
    5: 'too many exposures',
    6: 'multiple exposures with the same dither position',
    7: 'average seeing > maximum',
    8: 'exposures span more than one plugging.',
    9: 'set is incomplete but the plugging is no longer current.',
    10: 'status from database.'
}


def checkSet(set, flag=True, flagExposures=None, force=False, silent=False,
             **kwargs):
    """Checks if a set meets MaNGA's quality criteria.

    Checks the input Totoro.dbclasses.Set objects to confirm that it meets all
    the requirements set by MaNGA. If the set passes all checks, it is flagged
    in the DB with status `Good` or `Excellent` depending on its average
    seeing.

    Parameters
    ----------
    set : `Totoro.dbclasses.Set`
        The Totoro set object to be checked and flagged.

    flag : bool
        If True and `set` is not mock, the status of `set` will flagged `Good`
        or `Excellent` in the DB depending on the result of the check.

    flagExposures : bool or None
        The `flag` parameter to pass to `checkExposure` for each exposure in
        `set`. If None, the same value as in `flag` will be used.

    force : bool
        If True, the set will be reflagged even if a status has been previously
        assigned.

    Returns
    -------
    result : tuple
        A tuple in which the first element is the set status. This value can be
        one of `Incomplete`, `Good`, `Excellent`, `Bad`, `Unplugged`.
        The second element is the code error that caused the set to fail,
        if any.

    Error codes
    -----------

    0: no error.
    1: one or more exposures are invalid.
    2: HA range is greater than maximum allowed.
    3: seeing values out of range
    4: SN2 values out of range
    5: too many exposures
    6: multiple exposures with the same dither position
    7: average seeing > maximum
    8: exposures span more than one plugging.
    9: set is incomplete but the plugging is no longer current.
    10: status from database.

    """

    assert isinstance(set, Set), ('input set is not an instance of '
                                  'Totoro.dbclasses.Set')

    if set.isMock or any([exp.isMock for exp in set.totoroExposures]):
        flag = False
        pk = None
    else:
        pk = set.pk

    if flagExposures is None:
        flagExposures = flag

    # If the exposure is not mock, force is not True, and an override status
    # is set, returns that.
    if set.isMock is False and set.set_status_pk is not None:
        if set.status.label == 'Override Good':
            return ('Good', 10)
        elif set.status.label == 'Override Bad':
            return ('Bad', 10)
        elif ((set.status.label == 'Excellent' or
               set.status.label == 'Good') and not force):
            return ('Good', 10)
        elif set.status.label == 'Poor' and not force:
            return ('Bad', 10)
        else:
            pass

    # Empty set
    if len(set.totoroExposures) == 0:
        return ('Incomplete', 0)

    # Check if exposures are valid
    for exposure in set.totoroExposures:
        exposureCheck = exposure.checkExposure(flag=flagExposures, force=force)
        if exposureCheck[0] is False:
            message = ('set pk={0}: one or more exposures are invalid.'
                       .format(pk))
            return flagSet(set, 'Bad', 1, flag=flag, message=message,
                           silent=silent)

    # Checks range of observations
    HA = set.getHA()
    HALength = (HA[1] - HA[0]) % 360.
    maxHARange = config['set']['maxHARange']
    if HALength > maxHARange:
        message = ('set pk={0}: HA range is larger than {1} deg.'
                   .format(pk, maxHARange))
        return flagSet(set, 'Bad', 2, flag=flag, message=message,
                       silent=silent)

    # Checks seeing
    seeing = np.array([exp.seeing for exp in set.totoroExposures])
    maxSeeingRange = config['set']['maxSeeingRange']
    if np.max(seeing) - np.min(seeing) > maxSeeingRange:
        message = ('set pk={0} fails the seeing uniformity criteria'
                   .format(pk))
        return flagSet(set, 'Bad', 3, flag=flag, message=message,
                       silent=silent)

    # Checks SN2 uniformity
    sn2 = np.array([exp.getSN2Array() for exp in set.totoroExposures])
    maxSN2Factor = config['set']['maxSN2Factor']
    for ii in range(len(sn2)):
        for jj in range(ii, len(sn2)):
            sn2Ratio = sn2[ii] / sn2[jj]
            sn2Ratio[np.isnan(sn2Ratio)] = config['set']['maxSN2Factor']
            if (np.any(sn2Ratio > maxSN2Factor) or
                    np.any(sn2Ratio < (1. / maxSN2Factor))):
                message = ('set pk={0} fails the SN2 uniformity criteria'
                           .format(pk))
                return flagSet(set, 'Bad', 4, flag=flag, message=message,
                               silent=silent)

    # Checks dithers
    ditherPositions = config['set']['ditherPositions']
    setDitherPositions = np.array(set.getDitherPositions())

    if len(setDitherPositions) > len(ditherPositions):
        message = ('set pk={0} has {1} exposures!'
                   .format(pk, len(setDitherPositions)))
        return flagSet(set, 'Bad', 5, flag=flag, message=message,
                       silent=silent)

    # Checks if exposure dithers are unique
    if np.unique(setDitherPositions).size < setDitherPositions.size:
        message = ('set pk={0} has multiple exposures with '
                   'the same dither position'.format(pk))
        return flagSet(set, 'Bad', 6, flag=flag, message=message,
                       silent=silent)

    # Checks if all exposures belong to the same plugging.
    expPluggings = np.array([exp.getPlugging().pk
                             if exp.getPlugging() is not None else 0
                             for exp in set.totoroExposures])

    # Now compares the plugging pks.
    if len(expPluggings) > 0 and not np.all(expPluggings == expPluggings[0]):
        message = ('set pk={0} has exposures from different pluggings.'
                   .format(pk))
        return flagSet(set, 'Bad', 8, flag=flag, message=message,
                       silent=silent)

    # Checks if set is incomplete.
    if len(setDitherPositions) < len(ditherPositions):

        # Gets the current plugging
        activePlugging = None
        for exp in set.totoroExposures:
            if exp.getActivePlugging() is not None:
                activePlugging = exp.getActivePlugging().pk
                break

        allMock = all([exp.isMock for exp in set.totoroExposures])

        if not allMock and expPluggings[0] != activePlugging:
            message = 'set pk={0} is from an inactive plugging.'.format(pk)
            return flagSet(set, 'Unplugged', 9, flag=flag, message=message,
                           silent=silent)

        # Otherwise, the set is incomplete
        return flagSet(set, 'Incomplete', 0, flag=False, silent=silent)

    # Set is valid: assigns status
    goodSeeingLimit = config['set']['goodSeeing']
    # excellentSeeingLimit = config['set']['excellentSeeing']
    if np.mean(seeing) > goodSeeingLimit:
        message = ('set pk={0} has average seeing > {1:.1f}'
                   .format(pk, goodSeeingLimit))
        return flagSet(set, 'Bad', 7, flag=flag, message=message,
                       silent=silent)
    # elif np.mean(seeing) <= excellentSeeingLimit:
    #     return flagHelper('Excellent', 0)
    else:
        return flagSet(set, 'Good', 0, flag=flag, silent=silent)
