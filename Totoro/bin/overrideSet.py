#!/usr/bin/env python
# encoding: utf-8
"""
overrideSet.py

Created by José Sánchez-Gallego on 25 Sep 2015.
Licensed under a 3-clause BSD license.

Revision history:
    25 Sep 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division, print_function

import argparse
import os
import sys
import warnings
from builtins import map

import numpy as np
from sqlalchemy.orm.exc import NoResultFound

from Totoro import log
from Totoro.db import getConnection
from Totoro.dbclasses import fromPlateID
from Totoro.dbclasses.exposure import Exposure, setExposureStatus
from Totoro.dbclasses.plate_utils import getConsecutiveSets, removeOrphanedSets
from Totoro.dbclasses.set import Set, checkSet, setErrorCodes
from Totoro.exceptions import EmptySet, TotoroError, TotoroUserWarning


warnings.simplefilter('always')


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


def override(args):
    """Overrides a set of exposures as good or bad."""

    db = getConnection()
    session = db.Session()

    mode = args.mode.capitalize()
    exposures = args.EXPOSURE_NO
    verbose = args.verbose

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
        if verbose:
            log.info('creating new set pk={0}'.format(overridenSetPK))
        with session.begin():
            ss = db.mangaDB.Set(pk=overridenSetPK)
            session.add(ss)

    for totExp in totExposures:
        setExposureStatus(totExp, 'Override ' + mode)
        if verbose:
            log.info('changing exposure_no={0} mangaDB status to Override {1}'.format(
                totExp.exposure_no, mode))

        if totExp._mangaExposure.set_pk != overridenSetPK:
            with session.begin():
                totExp._mangaExposure.set_pk = overridenSetPK
            if verbose:
                log.info('changing set_pk for exposure_no={0} to {1}'.format(
                    totExp.exposure_no, overridenSetPK))

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

            if verbose:
                ss = Set(setPK, format='pk')
                newStatus = None if ss.status is None else ss.status.label
                if newStatus != origStatus:
                    log.info('set pk={0} status changed from {1} to {2}'.format(
                        setPK, origStatus, newStatus))

        except EmptySet:

            if verbose:
                log.info('set pk={0} is empty. Removing it.'.format(setPK))
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

    if verbose:
        log.info('changing status of set {0} to Override {1}'.format(overridenSetPK, mode))
        log.info('override was successful')

    removeOrphanedSets()

    return overridenSetPK


def removeStatus(args):
    """Removes set status."""

    db = getConnection()
    session = db.Session()

    setPK = args.SET_PK
    assert isinstance(setPK, int), 'SET_PK must be an integer'

    verbose = args.verbose
    reload = args.reload

    # Checks that the set exists
    with session.begin():
        ss = session.query(db.mangaDB.Set).get(setPK)
        if ss is None:
            raise TotoroError('set_pk={0} does not exist'.format(setPK))

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
            if verbose:
                log.info('removing set_pk and status for exposure_no {0}'
                         .format(exp.platedbExposure.exposure_no))
            exp.set_pk = None
            exp.exposure_status_pk = None

        # Removes the set
        session.delete(ss)

    # We run a removeOrphanedSets just to be sure
    if verbose:
        log.info('removing orphaned sets.')
    removeOrphanedSets()

    if reload:
        if verbose:
            log.info('reloading plateid {0}'.format(plateID))
        fromPlateID(plateID)

    warnings.warn('remember to check the set arrangement for plate_id={0} '
                  'after removing overridden sets.'.format(plateID))

    if verbose:
        log.info('set_pk={0} removed successfully'.format(setPK))

    return


def getInfo(args):
    """Returns information about a set composed by a list of exposures."""

    exposures = args.EXPOSURE_NO
    verbose = args.verbose

    if verbose:
        log.info('checking exposures ' + ', '.join(map(str, exposures)))

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


def main(argv=None):

    parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]))

    parser.add_argument(
        '--verbose', '-v', dest='verbose', action='store_true', help='Prints extra information.')

    subparsers = parser.add_subparsers(title='actions')

    parserGood = subparsers.add_parser(
        'good', help='overrides set as good.', description='Overrides set as good.')
    parserGood.add_argument(
        'EXPOSURE_NO',
        metavar='EXPOSURE_NO',
        type=int,
        nargs='*',
        help='The list of exposure_no(s) to be '
        'overridden as a new good set.')
    parserGood.set_defaults(func=override, mode='good')

    # parserBad = subparsers.add_parser('bad', help='overrides set as bad.',
    #                                   description='Overrides set as bad.')
    # parserBad.add_argument('EXPOSURE_NO', metavar='EXPOSURE_NO', type=int,
    #                        nargs='*', help='The list of exposure_no(s) to be '
    #                        'overridden as a new bad set.')
    # parserBad.set_defaults(func=override, mode='bad')

    parserInfo = subparsers.add_parser(
        'info',
        help='gets information about a set.',
        description='If the set exists it returns its real status. Otherwise, '
        'it returns status of the mock set and, if invalid, why it is so.')
    parserInfo.add_argument(
        'EXPOSURE_NO',
        metavar='EXPOSURE_NO',
        type=int,
        nargs='*',
        help='The list of exposure_no(s) to be tested.')
    parserInfo.set_defaults(func=getInfo)

    parserRemove = subparsers.add_parser(
        'remove', help='removes set status.', description='Removes set status.')
    parserRemove.add_argument(
        'SET_PK',
        metavar='SET_PK',
        type=int,
        help='The list of set_pk(s) of the sets for '
        'which the set status will be removed.')
    parserRemove.add_argument(
        '--reload',
        '-reload',
        action='store_true',
        help='If true, load the plate at which the set '
        'belongs after removing the set. '
        'This will force a new set arrangement of the '
        'exposures in the removed set.')
    parserRemove.set_defaults(func=removeStatus)

    if argv is not None:
        args = parser.parse_args(argv)
    else:
        args = parser.parse_args()

    return args.func(args)


if __name__ == '__main__':
    main()
