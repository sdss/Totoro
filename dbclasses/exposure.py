#!/usr/bin/env python
# encoding: utf-8
"""
exposure.py

Created by José Sánchez-Gallego on 13 Mar 2014.
Licensed under a 3-clause BSD license.

Revision history:
    13 Mar 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.internal.manga.Totoro.exceptions import TotoroError, NoMangaExposure
from sdss.internal.manga.Totoro import TotoroDBConnection
from sdss.internal.manga.Totoro import log, config, site, dustMap
from sdss.internal.manga.Totoro import utils
import numpy as np
from astropy import time
import warnings


db = TotoroDBConnection()
plateDB = db.plateDB
mangaDB = db.mangaDB
session = db.session


class Exposure(plateDB.Exposure):

    def __new__(cls, input=None, format='pk', parent='plateDB', *args,
                **kwargs):

        if input is None:
            return plateDB.Exposure.__new__(cls)

        base = cls.__bases__[0]

        if isinstance(input, base):
            instance = input

        elif isinstance(input, mangaDB.Exposure):
            instance = input.platedbExposure

        else:
            if parent.lower() == 'platedb':
                with session.begin():
                    instance = session.query(base).filter(
                        eval('{0}.{1} == {2}'.format(base.__name__,
                                                     format, input))).one()

            elif parent.lower() == 'mangadb':
                with session.begin():
                    instance = session.query(base).join(
                        mangaDB.Exposure).filter(
                            eval('mangaDB.Exposure.{0} == {1}'.format(
                                format, input))).one()

        instance.__class__ = cls

        return instance

    def __init__(self, input, format='pk', mock=False, *args, **kwargs):

        self._valid = None
        self._ditherPosition = None
        self._sn2Array = None
        self._seeing = None
        self._plugging = None
        self._mlhalimit = None
        self._haRange = None

        self.isMock = mock
        self.kwargs = kwargs

        self._mangaExposure = (self.mangadbExposure[0]
                               if len(self.mangadbExposure) > 0 else
                               mangaDB.Exposure())

        if self._mangaExposure.pk is None:
            warnings.warn('plateDB.Exposure.pk={0} has no mangaDB.Exposure '
                          'counterpart.'.format(self.pk), NoMangaExposure)

    def __repr__(self):
        return ('<Totoro Exposure (mangaDB.Exposure.pk={0}, exposure_no={1}, '
                'ditherPos={2}, valid={3})>'
                .format(self._mangaExposure.pk, self.exposure_no,
                        self.ditherPosition, self.valid))

    def update(self, **kwargs):
        """Reloads the exposure."""

        newSelf = Exposure(self.pk, fromat='pk', mock=self.isMock,
                           **self.kwargs)
        self = newSelf

        log.debug('Exposure pk={0} has been reloaded'.format(
                  self.pk))

    @classmethod
    def createMockExposure(cls, startTime=None, expTime=None,
                           ditherPosition=None, ra=None, dec=None,
                           plugging=None, silent=False, **kwargs):
        """Creates a mock exposure instance."""

        if ra is None or dec is None:
            raise TotoroError('ra and dec must be specified')

        newExposure = Exposure(None, mock=True, ra=ra, dec=dec, **kwargs)
        newExposure.pk = None if 'pk' not in kwargs else kwargs['pk']

        if startTime is None:
            startTime = time.Time.now().jd

        if expTime is None:
            expTime = config['exposure']['exposureTime']

        newExposure.ditherPosition = ditherPosition

        tt = time.Time(startTime, format='jd', scale='tai')
        t0 = time.Time(0, format='mjd', scale='tai')
        startTimePlateDB = (tt - t0).sec

        newExposure.start_time = startTimePlateDB
        newExposure.exposure_time = expTime
        newExposure._plugging = plugging

        newExposure.simulateObservedParamters(**kwargs)

        if not silent:
            log.debug('Created mock exposure with ra={0:.3f} and dec={0:.3f}'
                      .format(newExposure.ra, newExposure.dec))

        return newExposure

    def simulateObservedParamters(self, factor=None, dust=None, **kwargs):
        """Simulates the SN2 of the exposure, using dust extinction and airmass
        values."""

        self._seeing = 1.0

        if factor is None:
            factor = config['simulation']['factor']

        if dust is not None:
            self._dust = dust
        else:
            if dustMap is not None:
                self._dust = dustMap(self.ra, self.dec)
            else:
                self._dust = {'iIncrease': [1], 'gIncrease': [1]}

        haRange = self.getHA()
        ha = utils.calculateMean(haRange)
        self._airmass = utils.computeAirmass(self.dec, ha)

        nominalRed = config['simulation']['redSN2'] * factor
        nominalBlue = config['simulation']['blueSN2'] * factor
        alphaRed = config['simulation']['alphaRed']
        alphaBlue = config['simulation']['alphaBlue']

        sn2Red = (nominalRed / self._airmass ** alphaRed /
                  self._dust['iIncrease'][0])
        sn2Blue = (nominalBlue / self._airmass ** alphaBlue /
                   self._dust['gIncrease'][0])

        self._sn2Array = np.array([sn2Blue, sn2Blue, sn2Red, sn2Red])

        # self._valid = True
        # self.status = 'Good'

        return self._sn2Array

    @property
    def ra(self):
        """RA of the plate centre."""
        return self.getCoordinates()[0]

    @property
    def dec(self):
        """Dec of the plate centre."""
        return self.getCoordinates()[1]

    def getCoordinates(self):
        """Returns an array with the coordinates of the plate centre."""

        if 'ra' in self.kwargs and 'dec' in self.kwargs:
            if (self.kwargs['ra'] is not None and
                    self.kwargs['dec'] is not None):
                return np.array(
                    [self.kwargs['ra'], self.kwargs['dec']], np.float)
        else:
            pointing = (self.observation.plate_pointing.pointing)
            return np.array(
                [pointing.center_ra, pointing.center_dec], np.float)

    def getHA(self):
        """Returns the HA range in which the exposure was taken."""

        if self._haRange is not None:
            return self._haRange

        startTime = float(self.start_time)
        expTime = float(self.exposure_time)

        t0 = time.Time(0, format='mjd', scale='tai')
        tStart = t0 + time.TimeDelta(startTime, format='sec', scale='tai')

        lst = site.localSiderealTime(tStart.jd)
        ha0 = (lst * 15. - self.ra) % 360.

        ha = np.array([ha0, ha0 + expTime / 3600. * 15]) % 360.
        ha[ha > 180.] -= 360.

        self._haRange = ha

        return ha

    def getSN2Array(self, useNaN=True):
        """Returns an array with the SN2 of the exposure. The return
        format is [b1SN2, b2SN2, r1SN2, r2SN2]."""

        if self._sn2Array is not None:
            return self._sn2Array

        SN2Values = self._mangaExposure.sn2values[0]

        # None values are replaced with NaN, which are easier to work with
        # in numpy arrays.
        values = [SN2Values.b1_sn2, SN2Values.b2_sn2,
                  SN2Values.r1_sn2, SN2Values.r2_sn2]
        valuesNaN = [np.nan if value is None else value for value in values]

        if useNaN is True:
            return np.array(valuesNaN)
        else:
            return np.array(values)

    @property
    def valid(self):
        """Checks if an exposure is valid."""

        if self._valid is not None:
            return self._valid
        else:
            return self.isValid()[0]

    @valid.setter
    def valid(self, value):
        """Sets the validity of the exposure."""

        self._valid = value

    def isValid(self, flag=True, **kwargs):
        """Checks if an exposure is valid and the error code."""

        if self._valid is not None:
            return (self._valid, -1)

        status = checkExposure(self, flag=flag, **kwargs)

        if self.isMock:
            self._valid = status[0]

        return status

    def checkExposure(self, **kwargs):
        """Alias for `isValid`"""

        return self.isValid(**kwargs)

    @property
    def ditherPosition(self):
        """Gets the dither position for this exposure."""

        if self._ditherPosition is None and self.isMock is False:
            return self._mangaExposure.dither_position[0].upper()
        else:
            return self._ditherPosition

    @ditherPosition.setter
    def ditherPosition(self, value):
        """Sets the dither position for this exposure
        (does not write to the DB)"""

        self._ditherPosition = value

    @property
    def seeing(self):
        """Returns the seeing for this exposure."""

        if not self.isMock and self._seeing is None:
            return self._mangaExposure.seeing
        else:
            return self._seeing

    def getLST(self):
        """Returns the LST interval for this exposure."""

        ha0, ha1 = self.getHA()

        lst0 = (ha0 + self.ra) % 360. / 15
        lst1 = (ha1 + self.ra) % 360. / 15

        return np.array([lst0, lst1])

    def getUT(self, format=None):
        """Returns the UT interval in which this exposure was taken. If
        format is 'str', it returns a tuple with the UT0 and UT1 strings. If
        None or 'datetime', a tuple with datetime instances is returned."""

        startTime = float(self.start_time)
        t0 = time.Time(0, format='mjd', scale='tai')

        tStart = t0 + time.TimeDelta(startTime, format='sec', scale='tai')
        tEnd = tStart + time.TimeDelta(float(self.exposure_time), format='sec',
                                       scale='tai')

        ut0 = tStart.datetime
        ut1 = tEnd.datetime

        if format == 'str':
            return ('{0:%H:%M}'.format(ut0), '{0:%H:%M}'.format(ut1))
        else:
            return (ut0, ut1)

        return (ut0, ut1)

    def getJD(self):
        """Returns the JD interval in which this exposure was taken."""

        startTime = float(self.start_time)
        t0 = time.Time(0, format='mjd', scale='tai')

        tStart = t0 + time.TimeDelta(startTime, format='sec', scale='tai')
        tEnd = tStart + time.TimeDelta(
            float(self.exposure_time), format='sec', scale='tai')

        return (tStart.jd, tEnd.jd)

    def getPlatePK(self):
        """Returns the pk of the plate associated to this plate."""
        return int(self.observation.plate_pointing.plate.pk)

    def getPlateID(self):
        """Returns the plateid of the plate associated to this plate."""
        return int(self.observation.plate_pointing.plate.plate_id)

    def getMJD(self):
        """Gets the MJD for this exposure."""
        return self.mjd()

    def getPlugging(self):
        """Returns the Plugging instance from ModelClasses for the exposure."""

        if self._plugging is not None:
            return self._plugging

        if self.pk is None:
            return None
        else:
            # Caches the plugging information
            self._plugging = self.observation.plugging
            return self._plugging

    def getActivePlugging(self):
        """Gets the pk of the current plugging for a plate from one of its
        exposures."""

        if self.pk is None:
            return None

        plate = self.observation.plate_pointing.plate

        for plugging in plate.pluggings:
            if len(plugging.activePlugging) > 0:
                return plugging

        return None

    def isPlateDBValid(self):
        """Returns True if the exposure is Good or Override Good in plateDB."""

        if self.status.label in ['Good', 'Override Good']:
            return True
        return False

    @property
    def mlhalimit(self):
        """Returns the HA range for this exposure."""

        if self._mlhalimit is None:
            self._mlhalimit = utils.mlhalimit(self.dec)
        return self._mlhalimit


def flagExposure(exposure, status, errorCode, flag=True, message=None):
    """Helper function to log and flag exposures."""

    if message is not None:
        log.debug(message)

    if flag:
        if errorCode != 6:
            statusLabel = 'Totoro Good' if status else 'Totoro Bad'
            setExposureStatus(exposure, statusLabel)
        else:
            # If the exposure is not complete reduced, we remove the status,
            # if any.
            setExposureStatus(exposure, None)

    return (status, errorCode)


def setExposureStatus(exposure, status):
    """ Sets the status of an exposure.

    Parameters
    ----------
    exposure : `Totoro.dbclasses.Exposure`, `plateDB.Exposure`,
               `mangaDB.Exposure` or integer.
        The exposure that will receive the new status. Either a
        Totoro.Exposure, plateDB.Exposure or mangaDB.Exposure instance.
        Alternatively, the pk of the mangaDB.Exposure can be used.
    status : string or None
        The status to be set. It must be one of the values in
        mangaDB.ExposureStatus.label. If `status=None`, the status of the
        exposure is removed

    Returns
    -------
    result : bool
        Returns True if the status has been set correctly.

    Example
    -------
    To set the value of exposure pk=43 to "Totoro Good" ::
      >> setExposureStatus(43, 'Totoro Good')

    """

    if isinstance(exposure, Exposure):
        pk = exposure._mangaExposure.pk
    elif isinstance(exposure, db.plateDB.Exposure):
        pk = exposure.mangadbExposure[0].pk
    elif isinstance(exposure, db.mangaDB.Exposure):
        pk = exposure.pk
    else:
        pk = exposure

    # Gets the status_pk for the desired status and flags the exposure.
    with session.begin():

        try:
            if status is None:
                statusPK = None
            else:
                queryStatus = session.query(db.mangaDB.ExposureStatus).filter(
                    db.mangaDB.ExposureStatus.label == status).one()
                statusPK = queryStatus.pk
        except:
            raise TotoroError('status {0} not found in mangaDB.ExposureStatus'
                              .format(status))

        try:
            exp = session.query(db.mangaDB.Exposure).get(pk)
            exp.exposure_status_pk = statusPK
        except:
            raise TotoroError('failed while trying to set exposure_status for '
                              'mangaDB.Exposure={0}'.format(pk))

        # If status is Totoro Bad, removes set_pk
        if status == 'Totoro Bad':
            setPK = exp.set_pk
            if setPK is None:
                pass
            else:
                exp.set_pk = None
                ss = session.query(db.mangaDB.Set).get(setPK)
                if len(ss.exposures) == 0:
                    session.delete(ss)

    log.debug('mangaDB.Exposure.pk={0} set to {1}'.format(pk, status))

    return True


def checkExposure(exposure, flag=True, force=False, **kwargs):
    """Checks if a given exposures meets MaNGA's quality criteria.

    This function checks the validity of a MaNGA exposure and sets the
    appropriate status in mangaDB. If the exposure is mock, no flagging in the
    database is done, but the function still returns the status of the
    exposure. If the exposure has already been flagged with a status different
    than `Good`, the function returns that status, unless `force=True`.

    Parameters
    ----------
    exposure : `Totoro.dbclasses.Exposure`
        The Totoro exposure object to be checked and flagged.

    flag : bool
        If True and `exposure` is not mock, the status of `exposure` will
        be set to `Totoro Good` or `Totoro Bad` depending on the result of the
        check.

    force : bool
        If True, the exposure will be reflagged even if a status has been
        previously assigned.

    Returns
    -------
    result : tuple
        A tuple in which the first element is True if the exposure is valid and
        False otherwise; the second element is the code error that caused the
        exposure to fail.

    Error codes
    -----------

    0: no error.
    1: wrong dither position.
    2: exposure time too short.
    3: seeing too large.
    4: SN2 too low.
    5: HA range outside the range of visibility of the plate.
    6: exposure not completely reduced
    7: low transparency
    8: exposure taken during twilight
    9: plateDB status is Bad or Override Bad.
    10: status read from DB.
    """

    # Checks that exposure is a Totoro exposure.
    assert isinstance(exposure, Exposure), 'input is not a Totoro exposure.'

    pk = exposure._mangaExposure.pk

    # The following only applies to real exposure.
    if not exposure.isMock:
        # Gets the status of this exposure in plateDB and mangaDB
        if exposure._mangaExposure.status is None:
            mangaDBstatus = None
        else:
            mangaDBstatus = exposure._mangaExposure.status.label

        plateDBstatus = exposure.status.label

        # Now it does some checking to make sure these values are consistent.
        if plateDBstatus in ['Bad', 'Override Bad']:
            # Overrides Totoro status to bad, if needed
            if mangaDBstatus != 'Totoro Bad' or force:
                message = ('exposure is marked {0} in plateDB'
                           .format(plateDBstatus))
                return flagExposure(exposure, False, 10, flag=flag,
                                    message=message)
            else:
                return (False, 10)

        if mangaDBstatus is not None:
            if mangaDBstatus == 'Totoro Good' and not force:
                return (True, 10)
            elif mangaDBstatus == 'Totoro Bad' and not force:
                return (False, 10)
            elif mangaDBstatus == 'Override Good':
                return (True, 10)
            elif mangaDBstatus == 'Override Bad':
                return (False, 10)
            else:
                pass

    # If the exposure is mock, we set flag to False
    if exposure.isMock:
        flag = False

    # Checks plateDB status
    if not exposure.isMock and not exposure.isPlateDBValid():
        message = ('Invalid exposure. plateDB.Exposure.pk={0} '
                   'is marked {1} in plateDB.exposure_status.'
                   .format(pk, exposure.status.label))
        return flagExposure(exposure, False, 9, flag=flag, message=message)

    # Checks dither position
    validDitherPositions = config['exposure']['validDitherPositions']
    if exposure.ditherPosition not in validDitherPositions:
        message = ('Invalid exposure. plateDB.Exposure.pk={0} '
                   'has dither position {1}'.format(pk,
                                                    exposure.ditherPosition))
        return flagExposure(exposure, False, 1, flag=flag, message=message)

    # Checks transparency
    transparency = exposure._mangaExposure.transparency
    minTransparency = config['exposure']['transparency']
    if transparency is not None and transparency < minTransparency:
        message = ('Invalid exposure. plateDB.Exposure.pk={0} '
                   'has low transparency.'.format(pk))
        return flagExposure(exposure, False, 7, flag=flag, message=message)

    # Checks exposure time
    minExpTime = config['exposure']['minExpTime']
    expTime = exposure.exposure_time
    if expTime < minExpTime:
        message = ('Invalid exposure. plateDB.Exposure.pk={0} has an '
                   'exposure time shorter than the minimum acceptable.'
                   .format(pk))
        return flagExposure(exposure, False, 2, flag=flag, message=message)

    # Checks seeing
    maxSeeing = config['exposure']['maxSeeing']
    seeing = exposure.seeing
    if seeing > maxSeeing:

        message = ('Invalid exposure. plateDB.Exposure.pk={0} '
                   'has a seeing larger than the maximum acceptable.'
                   .format(pk))
        return flagExposure(exposure, False, 3, flag=flag, message=message)

    # Checks visibility window
    buffer = config['exposure']['exposureBuffer']
    exposureHA = exposure.getHA()
    visibilityWindow = np.array([-exposure.mlhalimit - buffer,
                                 exposure.mlhalimit + buffer])
    if not utils.isIntervalInsideOther(exposureHA, visibilityWindow,
                                       onlyOne=False, wrapAt=360):
        message = ('Invalid exposure. plateDB.Exposure.pk={0} '
                   'has HA range [{1}, {2}] that is outside the '
                   'visibility window of the plate [{3}, {4}]'
                   .format(pk, exposureHA[0], exposureHA[1],
                           visibilityWindow[0], visibilityWindow[1]))
        return flagExposure(exposure, False, 5, flag=flag, message=message)

    # Checks altitude of the object. Skipped if exposure is mock.
    if config['exposure']['checkTwilight'] is True and not exposure.isMock:
        # Avoids this check for mock exposures as it slows down simulations.
        maxSunAltitude = config['exposure']['maxSunAltitude']
        exposureDate = time.Time(exposure.getJD(), format='jd', scale='tai')
        sunAltitude = site.getSunAltitude(exposureDate)

        if np.any(sunAltitude > maxSunAltitude):
            message = ('Invalid exposure. plateDB.Exposure.pk={0}: '
                       'altitude the Sun at the time of the exposure '
                       '< {1:.1f} deg'.format(pk, maxSunAltitude))
            return flagExposure(exposure, False, 8, flag=flag, message=message)

    # Checks whether at least one camera for each arm has been reduced
    exposureSN2 = exposure.getSN2Array(useNaN=False)
    blue = exposureSN2[0:2]
    red = exposureSN2[2:]

    # If the exposure is partially reduced
    if not np.all(exposureSN2 >= 0):
        # If at least one camera in each spectrograph is reduced, does not flag
        # the exposure but returs True
        if np.any(blue >= 0) and np.any(red >= 0):
            message = ('plateDB.Exposure.pk={0}: not completely reduced but '
                       'temporarily considering it valid'.format(pk))
            return flagExposure(exposure, True, 6, flag=flag, message=message)
        else:
            # Otherwise, returns False
            message = ('plateDB.Exposure.pk={0}: not completely reduced but '
                       'temporarily considering it invalid'.format(pk))
            return flagExposure(exposure, False, 6, flag=flag, message=message)

    # Checks SN2
    minSN2red = config['SN2thresholds']['exposureRed']
    minSN2blue = config['SN2thresholds']['exposureBlue']

    blueSN2comp = blue < minSN2blue
    redSN2comp = red < minSN2red

    # If any of the cameras has a SN2 lower than the minimum, flag it as bad
    if any(blueSN2comp) or any(redSN2comp):
        message = ('Invalid exposure. plateDB.Exposure.pk={0} '
                   'has SN2 lower than the minimum acceptable.'.format(pk))
        return flagExposure(exposure, False, 4, flag=flag, message=message)

    # At this point, all the checks are done, so the exposure is valid.
    return flagExposure(exposure, True, 0, flag=flag)
