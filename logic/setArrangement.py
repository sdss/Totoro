#!/usr/bin/env python
# encoding: utf-8
"""
setArrangement.py

Created by José Sánchez-Gallego on 28 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    28 Aug 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import itertools
from mangaLogic import checkExposure
from sdss.internal.manga.Totoro.exceptions import NoMangaExposure
from sdss.internal.manga.Totoro import TotoroDBConnection, log, config, site
from sdss.internal.manga.Totoro.utils import intervals
from sqlalchemy import func
from scipy.misc import factorial
import collections
import warnings
import numpy as np
import time


def updatePlate(plate, **kwargs):
    """Finds missing exposures """

    log.debug('plate_id={0}: updating sets'.format(plate.plate_id))
    # print(plate)
    badSetStatus = checkBadSets(plate, **kwargs)
    if badSetStatus:
        plate = fixBadSets(plate)
        return True

    newExposures = getNewExposures(plate, **kwargs)

    if len(newExposures) > 0:

        log.info('plate_id={0}: {1} new exposure(s) found with '
                 'mangaDB.Exposure.pk={0}'.format(
                     plate.plate_id, len(newExposures),
                     ', '.join([str(exp.mangadbExposure[0].pk)
                                for exp in newExposures])))

        if len(newExposures) >= \
                config['setArrangement']['forceRearrangementMinExposures']:
            log.info('plate_id={0}: more than {1} new exposures found. '
                     'Triggering a complete set rearrangement.'
                     .format(plate.plate_id,
                             config['setArrangement']
                             ['forceRearrangementMinExposures']))
            updateStatus = rearrangeSets(plate, **kwargs)
        else:
            results = []
            for exp in newExposures:
                results.append(addExposure(exp, plate))
            if any(results):
                removeOrphanSets()
                updateStatus = True
            else:
                updateStatus = False
    else:
        log.debug('plate_id={0}: no new exposures found.'.format(
                  plate.plate_id))
        updateStatus = False

    if updateStatus:
        return True
    else:
        return False


def checkBadSets(plate, **kwargs):
    """Identifies bad sets and breaks them."""

    setQuality = [set.getQuality(silent=True)[0] for set in plate.sets]

    if 'Bad' in setQuality:
        nBadSets = len([True for ss in setQuality if ss == 'Bad'])
        log.debug('plate_id={0}: found {1} bad sets. Fixing them.'.format(
                  plate.plate_id, nBadSets))
        # plate = fixBadSets(plate, setQuality=setQuality)
        return True
    return False


def getValidSet(totoroExp, plate, setQuality=None):
    """Gets the best possible set for an exposure"""

    from sdss.internal.manga.Totoro import dbclasses

    if not totoroExp.isValid()[0]:
        return False

    if setQuality is None:
        setQuality = [set.getQuality(silent=True)[0] for set in plate.sets]

    incompleteSets = [set for nn, set in enumerate(plate.sets)
                      if setQuality[nn].lower() == 'incomplete']

    distanceToEndOfVisibilityWindow = [
        plate.getLSTRange()[1] - np.mean(set.getLSTRange())
        for set in incompleteSets]
    order = np.argsort(distanceToEndOfVisibilityWindow)
    incompleteSetsSorted = [incompleteSets[ii] for ii in order]

    mockSetQuality = []
    for set in incompleteSetsSorted:
        mockSet = dbclasses.Set.fromExposures(
            set.totoroExposures + [totoroExp], silent=True)
        mockSetQuality.append(mockSet.getQuality(silent=True)[0])
        if mockSetQuality[-1] in ['Good', 'Excellent']:
            return set

    try:
        firstIncompleteSet = mockSetQuality.index('Incomplete')
        return incompleteSetsSorted[firstIncompleteSet]
    except ValueError:
        return None


def addExposure(exp, plate):
    """Adds an exposure to a plate."""

    from sdss.internal.manga.Totoro import dbclasses

    if not isinstance(exp, dbclasses.Exposure):
        totoroExp = dbclasses.Exposure(exp, silent=True)
    else:
        totoroExp = exp

    validSet = getValidSet(exp, plate)
    if validSet is False:
        return False
    elif validSet is None:
        validSet = dbclasses.Set()
        plate.sets.append(validSet)

    validSet.totoroExposures.append(totoroExp)

    if plate.isMock:
        return True

    db = TotoroDBConnection()
    session = db.Session()

    with session.begin(subtransactions=True):
        if validSet.pk is None:
            validSet.pk = getNewSetPK()
            session.add(validSet)
            session.flush()
        totoroExp._mangaExposure.set_pk = validSet.pk
        session.add(totoroExp)
        session.flush()

    if validSet.pk is not None:
        log.info('adding mangaDB exposure pk={0} to set {1} -> {2} set.'
                 .format(totoroExp._mangaExposure.pk, validSet.pk,
                         validSet.getQuality(silent=True, flag=False)[0]))

    return True


def getNewSetPK():
    """Gets the lowest available set pk."""

    db = TotoroDBConnection()
    session = db.Session()

    with session.begin(subtransactions=True):
        setPKs = session.query(db.mangaDB.Set.pk).all()

    setPKs = np.unique(setPKs)
    return [xx for xx in np.arange(1, np.max(setPKs)+1) if xx not in setPKs][0]


def getNewExposures(plate, silent=True, **kwargs):
    """Gets unassigned exposures in a plate."""

    newExposures = []
    for exp in plate.getScienceExposures():
        if len(exp.mangadbExposure) == 0:
            warnings.warn('exposure pk={0} has no mangaDB '
                          'counterpart.'.format(exp.pk), NoMangaExposure)
            continue
        if checkExposure(exp, silent=silent, flag=True)[0] is False:
            continue
        if exp.mangadbExposure[0].set_pk is not None:
            continue
        newExposures.append(exp)

    return newExposures


def rearrangeSets(plate, **kwargs):
    """Assigns a set to each exposure for a given plate. Overwrites any current
    assignment"""

    log.info('plate_id={0}: rearranging sets'.format(plate.plate_id))

    t0 = time.time()
    optimalArrangement = getOptimalArrangement(plate, **kwargs)

    if optimalArrangement is False:
        return False

    for exposure in optimalArrangement['invalid']:
        removeSet(exposure._mangaExposure.set_pk)

    for set in optimalArrangement['sets']:
        updateSet(set)

    removeOrphanSets()

    t1 = time.time()
    log.info('plate_id={0}: rearrangement was successful (took {1:.5f}s).'
             .format(plate.plate_id, t1-t0))

    return True


def removeSet(set_pk, orphan=False):
    """Removes a set and sets set_pk to None for all its exposures."""

    if set_pk is None:
        return

    db = TotoroDBConnection()
    session = db.Session()

    with session.begin(subtransactions=True):
        set = session.query(db.mangaDB.Set).get(set_pk)
        if set is not None:
            if set.exposures is not None:
                for exp in set.exposures:
                    exp.set_pk = None
            session.delete(set)

    msg = 'removed orphan set pk={0}' if orphan else 'Removed set pk={0}'
    log.debug(msg.format(set_pk))

    return


def updateSet(set):
    """Updates the set_pk for the exposures in the set. Writes to the DB."""

    db = TotoroDBConnection()
    session = db.Session()

    setPKs = [exposure._mangaExposure.set_pk
              for exposure in set.totoroExposures]

    for ss in np.unique(setPKs):
        removeSet(ss)

    with session.begin(subtransactions=True):
        newSet = db.mangaDB.Set()
        newSet.pk = getNewSetPK()
        session.add(newSet)
        session.flush()
        for exposure in set.totoroExposures:
            mangaExp = exposure._mangaExposure
            mangaExp.set_pk = newSet.pk
            session.add(mangaExp)

    log.debug('created set pk={0} with mangaDB exposures pk={1}'
              .format(newSet.pk,
                      ', '.join([str(exp._mangaExposure.pk)
                                 for exp in set.totoroExposures])))
    return


def removeOrphanSets():
    """Removes sets without exposures."""

    db = TotoroDBConnection()
    session = db.Session()

    with session.begin(subtransactions=True):
        nullSets = session.query(db.mangaDB.Set).outerjoin(
            db.mangaDB.Exposure).group_by(db.mangaDB.Set.pk).having(
                func.count(db.mangaDB.Set.exposures) == 0).all()

    for nullSet in nullSets:
        removeSet(nullSet.pk, orphan=True)

    return


def getNumberPermutations(ditherPositions):
    """Estimates the number of permutations to check for a certain list of
    dithered positions."""

    if len(ditherPositions) == 0:
        return 0

    repDict = collections.defaultdict(int)
    for cc in ditherPositions:
        repDict[cc] += 1

    maxNDither = 0
    for key in repDict.keys():
        if repDict[key] > maxNDither:
            maxNDither = repDict[key]

    return int(factorial(maxNDither)**(len(repDict.keys())-1))


def getOptimalArrangement(plate, LST=None,
                          permutationLimit=None,
                          forceLimit=False, **kwargs):
    """Gets the best possible arrangement for the exposures in a plate."""

    silent = kwargs.pop('silent', True)

    from sdss.internal.manga.Totoro import dbclasses

    permutationLimit = config['setArrangement']['permutationLimit'] \
        if permutationLimit is None else permutationLimit

    exposures = [dbclasses.Exposure(exp.pk, parent='mangaDB', silent=True)
                 for exp in plate.getMangadbExposures()]

    validExposures = []
    invalidExposures = []
    for exposure in exposures:
        if checkExposure(exposure, forceReflag=True, flag=True,
                         silent=silent, **kwargs)[0]:
            validExposures.append(exposure)
        else:
            invalidExposures.append(exposure)

    if len(validExposures) == 0:
        log.info('no exposures to be rearranged')
        return False

    log.debug('rearranging sets with LST={0}'.format(LST))

    ditherPositions = [exp.ditherPosition for exp in validExposures]
    nPermutations = getNumberPermutations(ditherPositions)

    if nPermutations > permutationLimit:
        if forceLimit is False:
            log.info('plate_id={0}: hard limit for number of permutations '
                     'in rearrangement ({1}) reached. Not rearranging.'.format(
                         plate.plate_id, permutationLimit))
            return False
        else:
            log.info('plate_id={0}: hard limit for number of permutations '
                     'in rearrangement reached but ignoring because '
                     'forceLimit=True.'.format(plate.plate_id))

    permutations = calculatePermutations(ditherPositions)

    log.info('plate_id={0}: testing {1} combinations for {2} exposures. '
             'This might take a while.'
             .format(plate.plate_id, nPermutations, len(ditherPositions)))

    def getSetId(set):
        """Creates a unique identifier for a set based on the ids of its
        exposures."""
        return np.sum([id(exp) for exp in set.totoroExposures])

    plates = []
    setQuality = {}
    completion = []

    # Counts the actual number of permutations, from the iterator.
    # For testing purposes.
    permutationCounter = 0
    setRearrangementFactor = config['set']['setRearrangementFactor']

    for nn, permutation in enumerate(permutations):

        sets = []

        for setIndices in permutation:

            exposures = [validExposures[ii]
                         for ii in setIndices if ii is not None]

            set = dbclasses.Set.fromExposures(exposures, silent=True)
            sets.append(set)

            # To avoid calculating the state of a set more than one, creates
            # a dictionary with the quality of the set
            setId = getSetId(set)
            if setId not in setQuality:
                setQuality[setId] = set.getQuality(silent=True)[0]

        del set

        ra, dec = sets[-1].getCoordinates()
        mockPlate = dbclasses.Plate.fromSets(sets, ra=ra, dec=dec, dust=None,
                                             silent=True)

        # Instead of using Plate.getPlateCompletion, we calculate the plate
        # completion here using the setQuality dictionary. Way faster this
        # way.
        plateSN2 = np.nansum(
            [set.getSN2Array()
             if setQuality[getSetId(set)] in ['Excellent', 'Good']
             else np.array([0.0, 0.0, 0.0, 0.0])
             for set in sets], axis=0)

        blueSN2 = np.nanmean(plateSN2[0:2])
        blueCompletion = blueSN2 / config['SN2thresholds']['plateBlue']
        redSN2 = np.nanmean(plateSN2[2:])
        redCompletion = redSN2 / config['SN2thresholds']['plateRed']
        plateCompletion = np.min([blueCompletion, redCompletion])

        if (len(completion) == 0 or
                plateCompletion >=
                setRearrangementFactor * np.max(completion)):
            plates.append(mockPlate)
            completion.append(plateCompletion)

        del mockPlate

        if (nn+1)*100./nPermutations % 10 == 0:
            log.info('{0:d}% completed'.format(int((nn+1)*100./nPermutations)))

        permutationCounter += 1

    log.info('{0} permutations tested.'.format(permutationCounter))

    maxCompletion = np.max(completion)
    completion = np.array(completion)
    plates = np.array(plates)

    # Selects plates that have 0.9 the maximum completion or higher
    validPlates = plates[completion >= setRearrangementFactor * maxCompletion]
    validCompletion = completion[completion >=
                                 setRearrangementFactor * maxCompletion]
    sortCompletion = np.argsort(validCompletion)[::-1]
    del plates

    maxPlates = [validPlates[idx] for idx in sortCompletion]
    log.debug('{0} plates with completion > {1:.1f}*maxCompletion'.format(
              len(maxPlates), setRearrangementFactor))

    # For the selected plates, checks if they have bad sets and, in that case
    # breaks them into incomplete sets.
    for maxPlate in maxPlates:
        maxPlateSetQuality = [setQuality[getSetId(set)]
                              for set in maxPlate.sets]
        if 'Bad' in maxPlateSetQuality:
            maxPlate = fixBadSets(maxPlate, setQuality=maxPlateSetQuality)

    if len(maxPlates) == 1:
        optimumPlate = maxPlates[0]

    else:

        # If several plates have been selected in the previous step, selects
        # the one with the fewest number of incomplete sets.
        nIncompleteSets = []
        for tmpPlate in maxPlates:
            nn = 0
            for set in tmpPlate.sets:
                if set.getQuality(silent=True)[0] == 'Incomplete':
                    nn += 1
            nIncompleteSets.append(nn)

        minIncompleteSetPlates = [
            maxPlates[ii] for ii in
            np.where(np.array(nIncompleteSets) == np.min(nIncompleteSets))[0]]

        # If thre are no incomplete sets
        if np.min(nIncompleteSets) == 0:
            completion = [plateMin.getPlateCompletion()
                          for plateMin in minIncompleteSetPlates]
            optimumPlate = minIncompleteSetPlates[np.argmax(completion)]
        else:
            if len(minIncompleteSetPlates) == 0:
                optimumPlate = minIncompleteSetPlates[0]
            else:
                optimumPlate = getEarliestIncompletePlate(
                    minIncompleteSetPlates, LST=LST)

    return {'sets': optimumPlate.sets, 'invalid': invalidExposures}


def getEarliestIncompletePlate(incompletePlates, LST=None):
    """From a list of plates with incomplete sets, selects the one with the
    incomplete set with the earliest visibility window."""

    plateRA = incompletePlates[0].ra

    # Calculates the LST range of the plate and the current LST depending on
    # the value of the LST keyword.
    if LST is None or LST is True:
        plateLSTRange = incompletePlates[0].getLSTRange(intersect=True)
        LST = site.localSiderealTime()
    elif LST is False:
        plateLSTRange = incompletePlates[0].getLSTRange(intersect=False)
        LST = plateLSTRange[0]
    else:
        plateLSTRange = incompletePlates[0].getLSTRange(intersect=False)
        assert LST >= 0 and LST <= 24.  # Checks LST is a time.

    # If the desired LST is not withing the range of the plate, we use its
    # starting LST
    if not intervals.isPointInInterval(LST, plateLSTRange, wrapAt=24.):
        LST = plateLSTRange[0]

    HA = LST * 15. - plateRA
    if HA > 180:
        HA -= 360
    elif HA < -180:
        HA += 360

    minHAs = []
    for plate in incompletePlates:
        HAPlate = []
        for set in plate.sets:
            if set.getQuality(silent=True)[0] is not 'Incomplete':
                continue

            haMid = intervals.calculateMean(set.getHA(), wrapAt=360.)
            if haMid > 180:
                haMid -= 360.
            elif haMid < -180:
                haMid += 360
            HAPlate.append(haMid)
        minHAs.append(np.min(HAPlate))

    minHAs = np.array(minHAs)
    afterHA = minHAs[minHAs > HA]

    if len(afterHA) == 0:
        return incompletePlates[np.argmin(minHAs)]
    else:
        return incompletePlates[np.argmin(afterHA)]


def calculatePermutations(inputList):
    """Calculates all the possible permutations based on an input list of
    dithered positions."""

    pairs = [(nn, inputList[nn]) for nn in range(len(inputList))]
    pairs = sorted(pairs, key=lambda value: value[1])

    splitPairs = [list(bb) for aa, bb in itertools.groupby(
                  pairs, lambda value: value[1])]
    sortedPairs = sorted(splitPairs, key=lambda xx: len(xx))[::-1]

    indicesSeed = [[element[0] for element in sP] for sP in sortedPairs]
    nExpPerDither = [len(ii) for ii in indicesSeed]
    for ii in indicesSeed:
        if len(ii) < np.max(nExpPerDither):
            ii += [None] * (np.max(nExpPerDither)-len(ii))

    if len(indicesSeed) > 0:
        indices = [[tuple(indicesSeed[0])]]
        indices += [list(itertools.permutations(idx))
                    for idx in indicesSeed[1:]]
    else:
        indices = []

    cartesianProduct = itertools.product(*indices)
    for prod in cartesianProduct:
        yield list(itertools.izip_longest(*prod))


def fixBadSets(plate, setQuality=None):
    """Breaks bad sets into incomplete sets."""

    from sdss.internal.manga.Totoro import dbclasses

    if setQuality is None:
        setQuality = [set.getQuality(silent=True)[0] for set in plate.sets]

    exposuresToAssign = []
    setsToRemove = []
    newSetQuality = []

    for nn, set in enumerate(plate.sets):
        if setQuality[nn].lower() == 'bad':
            exposuresToAssign += set.totoroExposures
            setsToRemove.append(set)
        else:
            newSetQuality.append(setQuality[nn])

    map(plate.sets.remove, setsToRemove)

    for exposure in exposuresToAssign:
        validSet = getValidSet(exposure, plate, setQuality=newSetQuality)

        if validSet is False:
            continue

        elif validSet is None:
            validSet = dbclasses.Set.createMockSet(ra=plate.ra, dec=plate.dec,
                                                   silent=True)
            plate.sets.append(validSet)
            newSetQuality.append('Incomplete')

        validSet.totoroExposures.append(exposure)

    if not plate.isMock:
        map(removeSet, [set.pk for set in setsToRemove])
        map(updateSet, [set for set in plate.sets])

    return plate
