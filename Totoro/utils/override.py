#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Created: 2018-06-14
# @LastModified: 2018-06-14
# @Filename: override.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# @Copyright: José Sánchez-Gallego

import warnings

import numpy as np
from sqlalchemy.orm.exc import NoResultFound

from Totoro import log
from Totoro.db import getConnection
from Totoro.dbclasses import fromPlateID
from Totoro.dbclasses.exposure import Exposure, setExposureStatus
from Totoro.dbclasses.plate_utils import getConsecutiveSets, removeOrphanedSets
from Totoro.dbclasses.set import Set, checkSet, setErrorCodes
from Totoro.exceptions import EmptySet, TotoroError, TotoroUserWarning


def _getStatusPK(status):
    """Returns the pk for a set status label."""

    db = getConnection()
    session = db.Session()

    with session.begin():
        statuses = session.query(db.mangaDB.SetStatus).all()

    for st in statuses:
        if st.label.lower() == status.lower():
            return st.pk

    raise ValueError('set status {0} does not exist in the database'.format(status))


def _checkExposures(exposures):
    """Raises exceptions if exposures don't meed certain requirements."""

    if len(exposures) == 0:
        raise TotoroError('no exposures specified')
    elif len(exposures) > 3:
        raise TotoroError('sets must consist of <= 3 exposures')
    elif len(exposures) < 3:
        warnings.warn(
            'you are creating an overridden set with only '
            '{0} exposures'.format(len(exposures)), TotoroUserWarning)

    # Checks that exposures exist
    totExposures = []
    originalSetPKs = []
    plateIDs = []
    for expNo in exposures:
        try:
            totExp = Exposure(expNo, format='exposure_no', parent='plateDB')
        except NoResultFound:
            raise TotoroError('exposure_no {0} does not exist in the database'.format(expNo))

        if totExp._mangaExposure.pk is None:
            raise TotoroError('exposure_no {0} has no MaNGA counterpart.'.format(expNo))

        totExposures.append(totExp)
        originalSetPKs.append(totExp._mangaExposure.set_pk)
        plateIDs.append(totExp.getPlateID())

    return totExposures, originalSetPKs, plateIDs


def override(exposures, mode):
    """Overrides a set of exposures as good or bad."""

    db = getConnection()
    session = db.Session()

    overriddenStatusPK = _getStatusPK('Override ' + mode)

    totExposures, originalSetPKs, plateIDs = _checkExposures(exposures)

    if len(np.unique(plateIDs)) != 1:
        raise TotoroError('exposures belong to different plates.')
    else:
        plateID = plateIDs[0]

    # Stores plate completion before doing any change
    plate = fromPlateID(plateID, updateSets=False, fullCheck=False)
    preCompletion = plate.getPlateCompletion()

    # If the exposures are already a set, we flag it overridden good/bad.
    # If not, we create a new set with those exposures. We also change all
    # exposure statuses to overridden good/bad.
    if len(np.unique(originalSetPKs)) == 1 and originalSetPKs[0] is not None:
        overridenSetPK = originalSetPKs[0]
    else:
        overridenSetPK = getConsecutiveSets(1)[0]
        log.debug('creating new set pk={0}'.format(overridenSetPK))
        with session.begin():
            ss = db.mangaDB.Set(pk=overridenSetPK)
            session.add(ss)

    for totExp in totExposures:
        setExposureStatus(totExp, 'Override ' + mode)
        log.debug('changing exposure_no={0} mangaDB status to Override {1}'
                  .format(totExp.exposure_no, mode))

        if totExp._mangaExposure.set_pk != overridenSetPK:
            with session.begin():
                totExp._mangaExposure.set_pk = overridenSetPK
            log.debug('changing set_pk for exposure_no={0} to {1}'
                      .format(totExp.exposure_no, overridenSetPK))

    # Overrides the set good/bad.
    with session.begin():
        ss = session.query(db.mangaDB.Set).get(overridenSetPK)
        ss.set_status_pk = overriddenStatusPK

    # Reflags the remaining sets
    for setPK in np.unique(originalSetPKs):
        if setPK == overridenSetPK or setPK is None:
            continue

        try:
            ss = Set(setPK, format='pk')
            origStatus = None if ss.status is None else ss.status.label
            checkSet(ss, flag=True, force=True)

            ss = Set(setPK, format='pk')
            newStatus = None if ss.status is None else ss.status.label
            if newStatus != origStatus:
                log.debug('set pk={0} status changed from {1} to {2}'.format(
                    setPK, origStatus, newStatus))

        except EmptySet:

            log.debug('set pk={0} is empty. Removing it.'.format(setPK))
            with session.begin():
                ss = session.query(db.mangaDB.Set).get(setPK)
                session.delete(ss)
            continue

    # Checks if completion has changed for the plate and, if so, issues a
    # warning. Otherwise, issues a general warning.
    plate = fromPlateID(plateID, updateSets=False, fullCheck=False)
    postCompletion = plate.getPlateCompletion()

    if ((preCompletion < 1. and postCompletion > 1.) or
            (preCompletion > 1. and postCompletion < 1.)):
        warnings.warn(
            'plate completion has changed from {0:.2f} to {1:.2f}. '
            'Remember to check if the plate status is correct after '
            'overriding sets.'.format(preCompletion, postCompletion), TotoroUserWarning)
    else:
        warnings.warn('Remember to check if the plate status is correct after '
                      'overriding sets.', TotoroUserWarning)

    log.debug('changing status of set {0} to Override {1}'.format(overridenSetPK, mode))
    log.debug('override was successful')

    removeOrphanedSets()

    return overridenSetPK


def removeStatus(set_pk, reload=False):
    """Removes set status."""

    db = getConnection()
    session = db.Session()

    assert isinstance(set_pk, int), 'SET_PK must be an integer'

    # Checks that the set exists
    with session.begin():
        ss = session.query(db.mangaDB.Set).get(set_pk)
        if ss is None:
            raise TotoroError('set_pk={0} does not exist'.format(set_pk))

        # Checks that all exposures belong to the same plate
        plateIDs = []
        for exp in ss.exposures:
            totExp = Exposure(exp)
            plateIDs.append(totExp.getPlateID())

        if len(np.unique(plateIDs)) != 1:
            raise TotoroError('exposures belong to different plates.')
        else:
            plateID = plateIDs[0]

        # Removes the set assignment
        for exp in ss.exposures:
            log.debug('removing set_pk and status for exposure_no {0}'
                      .format(exp.platedbExposure.exposure_no))
            exp.set_pk = None
            exp.exposure_status_pk = None

        # Removes the set
        session.delete(ss)

    # We run a removeOrphanedSets just to be sure
    log.debug('removing orphaned sets.')
    removeOrphanedSets()

    if reload:
        log.debug('reloading plateid {0}'.format(plateID))
        fromPlateID(plateID)

    warnings.warn('remember to check the set arrangement for plate_id={0} '
                  'after removing overridden sets.'.format(plateID), TotoroUserWarning)

    log.debug('set_pk={0} removed successfully'.format(set_pk))

    return


def getInfo(exposures):
    """Returns information about a set composed by a list of exposures."""

    log.debug('checking exposures ' + ', '.join(map(str, exposures)))

    totExposures, originalSetPKs, plateIDs = _checkExposures(exposures)

    if len(np.unique(plateIDs)) != 1:
        raise TotoroError('exposures belong to different plates.')

    if len(np.unique(originalSetPKs)) == 1 and originalSetPKs[0] is not None:
        log.important('exposures belong to real set pk={0}'.format(originalSetPKs[0]))
        ss = Set(originalSetPKs[0])
    else:
        log.important('creating mock set for input exposures')
        ss = Set.fromExposures(totExposures)

    if ss.isMock is False and ss.set_status_pk is not None:
        status = ss.status.label
        code = 10
    else:
        status, code = ss.getStatus()

    log.important('set status is {0}'.format(status))
    log.important('error code is {0:d}: {1}'.format(code, setErrorCodes[code]))

    if code == 10:
        log.important('now ignoring status from database ... ')

        # If the set is overridden then all the exposures will also be
        # overridden. To make this test we change their status to None. This
        # won't be recorded in the DB.
        exposureStatusPK = []
        if 'Override' in status:
            for exp in totExposures:
                exposureStatusPK.append(exp._mangaExposure.exposure_status_pk)
                exp._mangaExposure.exposure_status_pk = None

        ss = Set.fromExposures(totExposures)
        statusMock, codeMock = ss.getStatus(flag=False, force=True, flagExposures=False)
        log.important('mock set status is {0}'.format(statusMock))
        log.important('error code is {0:d}: {1}'.format(codeMock, setErrorCodes[codeMock]))

        # Just in case, let's restore the exposure statuses
        if 'Override' in status:
            for ii, exp in enumerate(totExposures):
                exp._mangaExposure.exposure_status_pk = exposureStatusPK[ii]

    else:
        codeMock = statusMock = None

    if code not in [0, 9, 10] or codeMock not in [0, 9, 10]:
        warnings.warn(
            'this is not a comprehensive list of reasons why the '
            'set is invalid. Other conditions may be failing for '
            'this set apart from the specified here.', TotoroUserWarning)

    return (status, code, statusMock, codeMock)
