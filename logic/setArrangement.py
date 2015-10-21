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
from ..exceptions import TotoroError, NoMangaExposure
from .. import TotoroDBConnection
from .. import log, site, config
from sqlalchemy import func
from scipy.misc import factorial
import collections
import warnings
import numpy as np
from astropy.time import Time


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


def getValidSet(totoroExp, plate):

    from ..dbclasses import Set

    if not totoroExp.isValid()[0]:
        return False

    completeSets = []
    incompleteSets = []

    for set in plate.sets:
        quality = set.getQuality(flag=False)[0]
        if quality in ['Excellent', 'Good']:
            continue
        if quality == 'Bad':
            log.info('one bad set found (pk={0}). Triggering rearrangement.'
                     .format(set.pk))
            rearrangeSets(plate)  # This is probably not enough. Test.
            continue
        if quality == 'Incomplete':
            set.totoroExposures.append(totoroExp)
            tmpSetQuality = set.getQuality(flag=False)[0]
            if tmpSetQuality in ['Excellent', 'Good']:
                completeSets.append(set)
            elif tmpSetQuality == 'Incomplete':
                incompleteSets.append(set)
            set.totoroExposures.remove(totoroExp)

    newSet = Set.fromExposures([totoroExp], silent=True)

    if len(completeSets) > 0:
        validSet = completeSets[
            np.argmax([np.sum(set.getSN2Array()) for set in completeSets])]
    elif len(incompleteSets) > 0:
        validSet = incompleteSets[
            np.argmax([np.sum(set.getSN2Array()) for set in incompleteSets])]
    else:
        validSet = newSet

    return validSet


def addExposure(exp, plate):

    from ..dbclasses import Exposure

    if not isinstance(exp, Exposure):
        totoroExp = Exposure(exp, silent=True)
    else:
        totoroExp = exp

    validSet = getValidSet(exp, plate)
    if not validSet:
        return False

    db = TotoroDBConnection()
    session = db.Session()

    with session.begin():
        if validSet.pk is None:
            session.add(validSet)
            session.flush()
        totoroExp._mangaExposure.set_pk = validSet.pk
        session.add(totoroExp)

    if validSet.pk is not None:
        log.info('adding mangaDB exposure pk={0} to set {1} -> {2} set.'
                 .format(totoroExp._mangaExposure.pk, validSet.pk,
                         validSet.getQuality(flag=False)[0]))

    return True


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
        plate = Plate.fromSets(sets, ra=ra, dec=dec, dust=None, silent=True)
        plates.append(plate)
        del plate

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

    validPlates = plates[completion >= 0.9 * maxCompletion]
    validCompletion = completion[completion >= 0.9 * maxCompletion]
    sortCompletion = np.argsort(validCompletion)[::-1]

    maxPlates = [validPlates[idx] for idx in sortCompletion]
    log.debug('{0} plates with completion > 0.9*maxCompletion'.format(
              len(maxPlates)))

    if len(maxPlates) == 1:
        optimumPlate = maxPlates[0]

    else:

        incompleteMaxPlates = [createIncompleteSets(plate)
                               for plate in maxPlates]

        for plate in incompleteMaxPlates:
            for set in plate.sets:
                print(set.totoroExposures)
            print()

        # maxDitherIncomplete = [getMaxNDitherIncomplete(plate)
        #                        for plate in incompleteMaxPlates]

        # validIncomplete = [incompleteMaxPlates[ii]
        #                    for ii in np.where(maxDitherIncomplete ==
        #                                       np.max(maxDitherIncomplete))[0]]

        optimumPlate = getMinDistancePlate(validIncomplete)

    for set in optimumPlate.sets:
        print(set.totoroExposures)
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


def createIncompleteSets(plate):

    from ..dbclasses import Set

    while 'Bad' in [set.getQuality()[0] for set in plate.sets]:

        status = [set.getQuality()[0] for set in plate.sets]
        badSet = status.index('Bad')
        exposures = plate.sets[badSet].totoroExposures
        plate.sets.pop(badSet)

        for exp in exposures:
            newSet = getValidSet(exp, plate)
            # print(plate.sets[0].totoroExposures)
            print(exp, newSet.totoroExposures, id(newSet), [id(xx) for xx in plate.sets])
            if id(newSet) in [id(xx) for xx in plate.sets]:
                for set in plate.sets:
                    if id(set) == id(newSet):
                        plate.sets.remove(set)
                plate.sets.append(newSet)
            else:
                print('hi')
                plate.sets.append(newSet)


    # toRemove = []

    # for ii in range(len(plate.sets)):

    #     set = plate.sets[ii]

    #     if set.getQuality()[0] == 'Bad':

    #         toRemove.append(set)

    #         if len(set.totoroExposures) == 2:
    #             for exp in set.totoroExposures:
    #                 plate.sets.append(Set.fromExposures([exp], silent=True))

    #         elif len(set.totoroExposures) == 3:

    #             twoComb = False
    #             for expIter in itertools.combinations(set.totoroExposures, 2):

    #                 tmpSet = Set.fromExposures(expIter, silent=True)

    #                 if tmpSet.getQuality()[0] == 'Incomplete':

    #                     plate.sets.append(tmpSet)

    #                     for exp in set.totoroExposures:
    #                         if exp not in expIter:
    #                             plate.sets.append(
    #                                 Set.fromExposures([exp], silent=True))

    #                     twoComb = True

    #                     break

    #             if not twoComb:
    #                 for exp in set.totoroExposures:
    #                     plate.sets.append(Set.fromExposures([exp],
    #                                                         silent=True))

    #         else:
    #             raise TotoroError('Bad set found with either 1 or more than 3'
    #                               ' exposure. Something went wrong.')

    # for ss in toRemove:
    #     plate.sets.remove(ss)

    return plate


def getMaxNDitherIncomplete(plate):

    nDither = np.array([len(set.totoroExposures) for set in plate.sets])

    if np.min(nDither) == 3:
        return 3
    else:
        return np.max(nDither[nDither < 3])


def getMinDistancePlate(plates, startDate=None):

    distanceFromBeginningOfWindow = []

    for plate in plates:

        distances = []

        for set in plate.sets:

            if set.getQuality()[0] == 'Incomplete':
                distances.append(getDistance(set, startDate=startDate))

        minDistance = np.min(distances) if len(distances) > 0 else 0.
        distanceFromBeginningOfWindow.append(minDistance)

    optimumPlate = plates[np.argmin(distanceFromBeginningOfWindow)]

    return optimumPlate


def getDistance(set, startDate):

    if startDate is None:
        startDate = Time.now().jd

    nExpNeeded = 3 - len(set.totoroExposures)
    minLength = nExpNeeded * config['exposure']['exposureTime'] / 3600.

    setLST = set.getLSTRange()
    startDateLST = site.localSiderialTime(startDate)

    distance = setLST[1] - startDateLST

    if distance < 0:
        return distance % 24
    else:
        if distance > minLength:
            return distance % 24
        else:
            return (setLST[0] - startDateLST) % 24
