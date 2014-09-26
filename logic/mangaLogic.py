#!/usr/bin/env python
# encoding: utf-8
"""
mangaLogic.py

Created by José Sánchez-Gallego on 4 Jul 2014.
Licensed under a 3-clause BSD license.

Contains functions and classes related with the scheduling logic for
MaNGA plates.

Revision history:
    4 Jul 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.internal.manga.Totoro import TotoroDBConnection, config, log
import numpy as np
from sdss.internal.manga.Totoro import utils


def removeSet(set_pk):
    """Removes a set."""

    db = TotoroDBConnection()
    session = db.Session()

    with session.begin(subtransactions=True):
        set = session.query(db.mangaDB.Set).get(set_pk)
        if set is None:
            log.debug('removing set pk={0} failed because the set does not'
                      ' exist.'.format(set_pk))
            return False
        else:
            session.delete(set)
            log.debug('removed set pk={0}'.format(set_pk))
            return True


def checkExposure(exposure, format='pk', parent='plateDB', flag=True,
                  silent=False, forceReflag=False, **kwargs):
    """Checks if a given exposures meets MaNGA's quality criteria.

    Error codes:

    0: no error.
    1: wrong dither position.
    2: exposure time too short.
    3: seeing too large.
    4: SN2 too low.
    5: HA range outside the range of visibility of the plate.
    6: exposure not completely reduced
    10: manually overriden bad.
    """

    from sdss.internal.manga.Totoro import dbclasses

    def flagHelper(status, errorCode, message=None):
        """Helper function to log and flag exposures."""
        if not silent and message is not None:
            log.debug(message)
        if flag:
            statusLabel = 'Totoro Good' if status else 'Totoro Bad'
            setExposureStatus(exposure, statusLabel)
        return (status, errorCode)

    # Performs a quick check to see if the exposure is alredy flagged.
    statusLabel = None
    if isinstance(exposure, dbclasses.Exposure):
        if not exposure.isMock and exposure._mangaExposure.status is not None:
            statusLabel = exposure._mangaExposure.status.label
        elif exposure.isMock:
            flag = False
    else:
        if hasattr(exposure, 'mangadbExposure'):
            if exposure.mangadbExposure[0].status is not None:
                statusLabel = exposure.mangadbExposure[0].status.label
        else:
            exposure = dbclasses.Exposure(exposure, format=format,
                                          parent=parent, silent=True)
            statusLabel = exposure._mangaExposure.status.label

    if statusLabel is not None:
        # If the exposure has been overriden good or bad outside Totoro,
        # we don't want to change that.
        if statusLabel.lower() == 'override good':
            return (True, 10)
        elif statusLabel.lower() in ['override bad', 'bad']:
            return (False, 10)
        elif not forceReflag and statusLabel.lower() == 'totoro good':
            return (True, 10)
        elif not forceReflag and statusLabel.lower() == 'totoro bad':
            return (False, 10)

    # If the exposure is not flagged, perform the QA tests
    if not isinstance(exposure, dbclasses.Exposure):
        exposure = dbclasses.Exposure(exposure, format=format,
                                      parent=parent, silent=silent)

    # Checks if the exposure has been completely reduce. If not, returns False
    # but does not flag it as Totoro Good/Bad.
    if None in list(exposure.getSN2Array()):
        setExposureStatus(exposure, 'Good')
        return (False, 6)

    # Checks dither position
    if (exposure.ditherPosition not in
            config['exposure']['validDitherPositions']
            and exposure.ditherPosition is not None):
        return flagHelper(False, 1,
                          'Invalid exposure. plateDB.Exposure.pk={0} '
                          'has dither position {1}'
                          .format(exposure._mangaExposure.pk,
                                  exposure.ditherPosition))

    # Checks exposure time
    minExpTime = config['exposure']['minExpTime']
    expTime = exposure.exposure_time
    if expTime < minExpTime:
        return flagHelper(False, 2,
                          'Invalid exposure. plateDB.Exposure.pk={0} has an '
                          'exposure time shorter than the minimum acceptable.'
                          .format(exposure._mangaExposure.pk))

    # Checks seeing
    maxSeeing = config['exposure']['maxSeeing']
    seeing = exposure.seeing
    if seeing > maxSeeing:
        return flagHelper(False, 3,
                          'Invalid exposure. plateDB.Exposure.pk={0} '
                          'has a seeing larger than the maximum acceptable.'
                          .format(exposure._mangaExposure.pk))

    # Checks SN2
    minSN2red = config['SN2thresholds']['exposureRed']
    minSN2blue = config['SN2thresholds']['exposureBlue']
    snArray = exposure.getSN2Array()

    if np.any(snArray[0:2] < minSN2blue) or np.any(snArray[2:] < minSN2red):
        return flagHelper(False, 4,
                          'Invalid exposure. plateDB.Exposure.pk={0} has SN2 '
                          'lower than the minimum acceptable.'
                          .format(exposure._mangaExposure.pk))

    # Checks visibility window
    visibilityWindow = np.array([-exposure.mlhalimit, exposure.mlhalimit])
    HA = exposure.getHA()
    if not utils.isIntervalInsideOther(HA, visibilityWindow,
                                       onlyOne=False, wrapAt=360):
        return flagHelper(False, 5,
                          'Invalid exposure. plateDB.Exposure.pk={0} '
                          'has HA range [{1}, {2}] that is outside the '
                          'visibility window of the plate [{3}, {4}]'
                          .format(exposure._mangaExposure.pk, HA[0], HA[1],
                                  visibilityWindow[0], visibilityWindow[1]))

    # If we are here is that everything went ok
    return flagHelper(True, 0)

    return


def setExposureStatus(exposure, status, **kwargs):
    """ Sets the status of an exposure.

    Parameters
    ----------
    exposure : Totoro.Exposure, plateDB.Exposure, mangaDB.Exposure or int
        The exposure that will receive the new status. Either a
        Totoro.Exposure, plateDB.Exposure or mangaDB.Exposure instance.
        Alternatively, the pk of the mangaDB.Exposure can be used.
    status : string
        The status to be set. It must be one of the values in
        mangaDB.ExposureStatus.label.

    Returns
    -------
    result : bool
        Returns True if the status has been set correctly.

    Example
    -------
    To set the value of exposure pk=43 to "Totoro Good" ::
      >> setExposureStatus(43, 'Totoro Good')

    """

    from sdss.internal.manga.Totoro import dbclasses

    db = TotoroDBConnection()
    session = db.Session()

    if isinstance(exposure, dbclasses.Exposure):
        pk = exposure._mangaExposure.pk
    elif isinstance(exposure, db.plateDB.Exposure):
        pk = exposure.mangadbExposure[0].pk
    elif isinstance(exposure, db.mangaDB.Exposure):
        pk = exposure.pk
    else:
        pk = exposure

    with session.begin(subtransactions=True):
        queryStatus = session.query(db.mangaDB.ExposureStatus).filter(
            db.mangaDB.ExposureStatus.label == status).one()
        statusPK = queryStatus.pk
        exp = session.query(db.mangaDB.Exposure).get(pk)
        exp.exposure_status_pk = statusPK

    log.debug('mangaDB.Exposure.pk={0} set to {1}'.format(exposure.pk, status))

    return True


def checkSet(input, flag=True, flagExposures=True, silent=False,
             midPoint=None, forceReflag=False, **kwargs):
    """Checks if a set meets MaNGA's quality criteria. Returns one of the
    following values: 'Good', 'Excellent', 'Poor', 'Bad', 'Incomplete'.

    Error codes:

    0: no error.
    1: one or more exposures are invalid.
    2: HA range is greater than maximum allowed.
    3: seeing values out of range
    4: SN2 values out of range
    5: too many exposures
    6: multiple exposures with the same dither position
    7: average seeing > maximum
    10: from set status.

    """

    from sdss.internal.manga.Totoro import dbclasses

    if isinstance(input, dbclasses.Set):
        set = input
        if set.isMock:
            flag = False
    else:
        set = dbclasses.Set(input, silent=silent)

    if set.isMock is False:
        if set.set_status_pk is not None and not forceReflag:
            return (set.status.label, 10)

    def flagHelper(statusLabel, errorCode, message=None):
        """Helper function to log and flag sets."""
        if not silent and message is not None:
            log.debug(message)
        if flag:
            setSetStatus(set, statusLabel)

        return (statusLabel, errorCode)

    if len(set.totoroExposures) == 0:
        return ('Incomplete', 0)

    # Check if exposures are valid
    for exposure in set.totoroExposures:
        exposureCheck = checkExposure(exposure, flag=flagExposures)
        if exposureCheck[0] is False:
            return flagHelper('Bad', 1,
                              'set pk={0}: one or more exposures are invalid.'
                              .format(set.pk))

    # Checks range of observations
    HA = set.getHA(midPoint=midPoint)
    HALength = (HA[1] - HA[0]) % 360.
    if HALength > config['set']['maxHARange']:
        return flagHelper('Bad', 2,
                          'set pk={0}: HA range is larger than {1} deg.'
                          .format(set.pk, config['set']['maxHARange']))

    # Checks seeing
    seeing = np.array([exp.seeing for exp in set.totoroExposures])
    if np.max(seeing) - np.min(seeing) > config['set']['maxSeeingRange']:
        return flagHelper('Bad', 3,
                          'set pk={0} '.format(set.pk) +
                          'fails the seeing uniformity criteria')

    # Checks SN2 uniformity
    sn2 = np.array([exp.getSN2Array() for exp in set.totoroExposures])
    for ii in range(len(sn2)):
        for jj in range(ii, len(sn2)):
            sn2Ratio = sn2[ii] / sn2[jj]
            if np.any(sn2Ratio > config['set']['maxSN2Factor']) or \
                    np.any(sn2Ratio < (1. / config['set']['maxSN2Factor'])):
                return flagHelper('Bad', 4,
                                  'set pk={0} '.format(set.pk) +
                                  'fails the SN2 uniformity criteria')

    # Checks dithers
    ditherPositions = config['set']['ditherPositions']
    setDitherPositions = np.array(set.getDitherPositions())

    if len(setDitherPositions) > len(ditherPositions):
        return flagHelper('Bad', 5,
                          'set pk={0} has {1} exposures!'
                          .format(set.pk, len(setDitherPositions)))

    if not all([sD is None for sD in setDitherPositions]):
        if np.unique(setDitherPositions).size < setDitherPositions.size:
            return flagHelper('Bad', 6,
                              'set pk={0} has multiple exposures with '
                              'the same dither position'.format(set.pk))

    # Checks if set is incomplete
    if len(setDitherPositions) < len(ditherPositions):
        return flagHelper('Incomplete', 0,
                          'set pk={0} is incomplete.'.format(set.pk))

    # Set is valid: assigns status
    if np.mean(seeing) > config['set']['goodSeeing']:
        return flagHelper('Bad', 7)
    elif np.mean(seeing) <= config['set']['excellentSeeing']:
        return flagHelper('Excellent', 0)
    else:
        return flagHelper('Good', 0)


def setSetStatus(set, status):
    """Sets the status of a set."""

    from sdss.internal.manga.Totoro import dbclasses

    db = TotoroDBConnection()
    session = db.Session()

    if isinstance(set, (dbclasses.Set, db.mangaDB.Set)):
        pk = set.pk
    else:
        pk = set

    with session.begin(subtransactions=True):
        try:
            queryStatus = session.query(db.mangaDB.SetStatus).filter(
                db.mangaDB.SetStatus.label == status).one()
            statusPK = queryStatus.pk
        except:
            statusPK = None

        ss = session.query(db.mangaDB.Set).get(pk)
        ss.set_status_pk = statusPK

    log.debug('mangaDB.Set.pk={0} set to {1}'.format(pk, status))

    return True
