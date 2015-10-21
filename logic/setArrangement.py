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
from .mangaLogic import checkExposure
from ..exceptions import NoMangaExposure
from .. import TotoroDBConnection
from .. import log, config
from sqlalchemy import func
from scipy.misc import factorial
import collections
import warnings
import numpy as np


def updateSets(plate, **kwargs):
    """Finds missing exposures """

    log.debug('updating sets for plate_id={0}'.format(plate.plate_id))

    newExposures = getNewExposures(plate)

    if len(newExposures) > 0:

        log.info('{0} new exposure(s) found with '
                 'mangaDB.Exposure.pk={1}'.format(
                     len(newExposures),
                     ', '.join([str(exp.mangadbExposure[0].pk)
                                for exp in newExposures])))

        if len(newExposures) >= \
                config['setArrangement']['forceRearrangementMinExposures']:
            log.info('more than {0} new exposures found. Triggering a complete'
                     ' set rearrangement.'
                     .format(config['setArrangement']
                             ['forceRearrangementMinExposures']))
            return rearrangeSets(plate, **kwargs)
        else:
            results = []
            for exp in newExposures:
                results.append(addExposure(exp, plate))
            if any(results):
                removeOrphanSets()
                return True
            else:
                return False
    else:
        log.debug('no new exposures found.')
        return False


def getValidSet(totoroExp, plate, setStatus=None):

    from ..dbclasses import Set

    if not totoroExp.isValid()[0]:
        return False

    if setStatus is None:
        setStatus = [set.getQuality()[0] for set in plate.sets]

    incompleteSets = [set for nn, set in enumerate(plate.sets)
                      if setStatus[nn].lower() == 'incomplete']

    distanceToEndOfVisibilityWindow = [
        plate.getLSTRange()[1] - np.mean(set.getLSTRange())
        for set in incompleteSets]
    order = np.argsort(distanceToEndOfVisibilityWindow)
    incompleteSetsSorted = [incompleteSets[ii] for ii in order]

    mockSetQuality = []
    for set in incompleteSetsSorted:
        mockSet = Set.fromExposures(set.totoroExposures + [totoroExp],
                                    silent=True)
        mockSetQuality.append(mockSet.getQuality()[0])
        if mockSetQuality[-1] in ['Good', 'Excellent']:
            return set

    try:
        firstIncompleteSet = mockSetQuality.index('Incomplete')
        return incompleteSetsSorted[firstIncompleteSet]
    except ValueError:
        return None


def addExposure(exp, plate):

    from ..dbclasses import Exposure, Set

    if not isinstance(exp, Exposure):
        totoroExp = Exposure(exp, silent=True)
    else:
        totoroExp = exp

    validSet = getValidSet(exp, plate)
    if validSet is False:
        return False
    elif validSet is None:
        validSet = Set()
        plate.sets.append(validSet)

    validSet.totoroExposures.append(totoroExp)

    db = TotoroDBConnection()
    session = db.Session()

    with session.begin():
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
                         validSet.getQuality(flag=False)[0]))

    return True


def getNewSetPK():
    """Gets the lowest available set pk."""

    db = TotoroDBConnection()
    session = db.Session()

    with session.begin(subtransactions=True):
        setPKs = session.query(db.mangaDB.Set.pk).all()

    setPKs = np.unique(setPKs)
    return [xx for xx in np.arange(1, np.max(setPKs)+1) if xx not in setPKs][0]


def getNewExposures(plate):

    newExposures = []
    for exp in plate.getScienceExposures():
        if len(exp.mangadbExposure) == 0:
            warnings.warn('exposure pk={0} has no mangaDB '
                          'counterpart.'.format(exp.pk), NoMangaExposure)
            continue
        if exp.mangadbExposure[0].set_pk is not None:
            continue
        if checkExposure(exp, silent=True, flag=True)[0] is False:
            continue
        newExposures.append(exp)

    return newExposures


def rearrangeSets(plate, **kwargs):
    """Assigns a set to each exposure for a given plate. Overwrites any current
    assignment"""

    log.info('reorganising sets for plate_id={0}'.format(plate.plate_id))

    optimalArrangement = getOptimalArrangement(plate, **kwargs)

    if optimalArrangement is False:
        return False

    for exposure in optimalArrangement['invalid']:
        removeSet(exposure._mangaExposure.set_pk)

    for set in optimalArrangement['sets']:
        updateSet(set)

    removeOrphanSets()

    log.info('rearrangement was successful')

    return True


def removeSet(set_pk, orphan=False):

    if set_pk is None:
        return

    db = TotoroDBConnection()
    session = db.Session()

    with session.begin():
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

    db = TotoroDBConnection()
    session = db.Session()

    setPKs = [exposure._mangaExposure.set_pk
              for exposure in set.totoroExposures]

    for ss in setPKs:
        removeSet(ss)

    with session.begin():
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

    db = TotoroDBConnection()
    session = db.Session()

    with session.begin():
        nullSets = session.query(db.mangaDB.Set).outerjoin(
            db.mangaDB.Exposure).group_by(db.mangaDB.Set.pk).having(
                func.count(db.mangaDB.Set.exposures) == 0).all()

    for nullSet in nullSets:
        removeSet(nullSet.pk, orphan=True)

    return


def getNumberPermutations(ditherPositions):

    repDict = collections.defaultdict(int)
    for cc in ditherPositions:
        repDict[cc] += 1

    return int(
        np.product(
            [factorial(ii) for ii in sorted(repDict.values())[1:]]))


def getOptimalArrangement(plate, startDate=None,
                          expLimit=config['setArrangement']['exposureLimit'],
                          forceLimit=False, **kwargs):

    from ..dbclasses import Exposure, Set, Plate

    exposures = [Exposure(exp.pk, parent='mangaDB', silent=True)
                 for exp in plate.getMangadbExposures()]

    validExposures = []
    invalidExposures = []
    for exposure in exposures:
        if checkExposure(exposure, **kwargs)[0]:
            validExposures.append(exposure)
        else:
            invalidExposures.append(exposure)

    if len(validExposures) == 0:
        return False

    if len(validExposures) > expLimit:
        if forceLimit is False:
            log.info('hard limit for number of exposures in rearrangement '
                     '({0}) reached. Not rearranging.'.format(expLimit))
            return False
        else:
            log.debug('hard limit for number of exposures in rearrangement '
                      'reached but ignoring because forceLimit=True.')

    ditherPositions = [exp.ditherPosition for exp in validExposures]
    permutations = calculatePermutations(ditherPositions)

    log.info('testing {0} combinations. This might take a while.'.format(
        getNumberPermutations(ditherPositions)))

    def getSetId(set):
        """Creates a unique identifier for a set based on the ids of its
        exposures."""
        return np.sum([id(exp) for exp in set.totoroExposures])

    plates = []
    setQuality = {}
    completion = []
    for permutation in permutations:

        sets = []

        for setIndices in permutation:

            exposures = [validExposures[ii]
                         for ii in setIndices if ii is not None]
            set = Set.fromExposures(exposures, silent=True)
            sets.append(set)

            # To avoid calculating the state of a set more than one, creates
            # a dictionary with the quality of the set
            setId = getSetId(set)
            if setId not in setQuality:
                setQuality[setId] = set.getQuality()[0]

        del set

        ra, dec = sets[-1].getCoordinates()
        mockPlate = Plate.fromSets(sets, ra=ra, dec=dec, dust=None,
                                   silent=True)

        plates.append(mockPlate)
        del mockPlate

        # Instead of using Plate.getPlateCompletion, we calculate the plate
        # completion here using the setQuality dictionary. Way faster this
        # way.
        plateSN2 = np.sum(
            [set.getSN2Array()
             if setQuality[getSetId(set)] in ['Excellent', 'Good']
             else np.array([0.0, 0.0, 0.0, 0.0])
             for set in sets], axis=0)
        plateSN2[0:2] /= config['SN2thresholds']['plateBlue']
        plateSN2[2:] /= config['SN2thresholds']['plateRed']

        completion.append(np.min(plateSN2))

    maxCompletion = np.max(completion)
    completion = np.array(completion)
    plates = np.array(plates)

    # Selects plates that have 0.9 the maximum completion or higher
    validPlates = plates[completion >= 0.9 * maxCompletion]
    validCompletion = completion[completion >= 0.9 * maxCompletion]
    sortCompletion = np.argsort(validCompletion)[::-1]

    maxPlates = [validPlates[idx] for idx in sortCompletion]
    log.debug('{0} plates with completion > 0.9*maxCompletion'.format(
              len(maxPlates)))

    # For the selected plates, checks if they have bad sets and, in that case
    # breaks them into incomplete sets.
    for maxPlate in maxPlates:
        setStatus = [setQuality[getSetId(set)] for set in maxPlate.sets]
        if 'Bad' in setStatus:
            fixBadSets(maxPlate, setStatus=setStatus)

    if len(maxPlates) == 1:
        optimumPlate = maxPlates[0]

    else:

        # If several plates have been selected in the previous step, selects
        # the one with the fewest number of incomplete sets.
        nIncompleteSets = []
        for tmpPlate in maxPlates:
            nn = 0
            for set in tmpPlate.sets:
                if set.getQuality()[0] == 'Incomplete':
                    nn += 1
            nIncompleteSets.append(nn)

        optimumPlate = maxPlates[np.argmin(nIncompleteSets)]

    return {'sets': optimumPlate.sets, 'invalid': invalidExposures}


def calculatePermutations(inputList):

    pairs = [(nn, inputList[nn]) for nn in range(len(inputList))]
    pairs = sorted(pairs, key=lambda value: value[1])

    splitPairs = [list(bb) for aa, bb in itertools.groupby(
                  pairs, lambda value: value[1])]

    indices = [[element[0] for element in sP] for sP in splitPairs]

    indices[0] = [indices[0]]
    indices[1:] = [list(itertools.permutations(idx)) for idx in indices[1:]]

    cartesianProduct = itertools.product(*indices)
    for prod in cartesianProduct:
        yield list(itertools.izip_longest(*prod))


def fixBadSets(plate, setStatus=None):

    from ..dbclasses import Set

    if setStatus is None:
        setStatus = [set.getQuality()[0] for set in plate.sets]

    exposuresToAssign = []
    setsToRemove = []
    newSetStatus = []

    for nn, set in enumerate(plate.sets):
        if setStatus[nn].lower() == 'bad':
            exposuresToAssign += set.totoroExposures
            setsToRemove.append(set)
        else:
            newSetStatus.append(setStatus[nn])

    map(plate.sets.remove, setsToRemove)

    for exposure in exposuresToAssign:
        validSet = getValidSet(exposure, plate, setStatus=newSetStatus)

        if validSet is None:
            validSet = Set.createMockSet(ra=plate.ra, dec=plate.dec,
                                         silent=True)
            plate.sets.append(validSet)
            newSetStatus.append('Incomplete')

        validSet.totoroExposures.append(exposure)

    return plate
