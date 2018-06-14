#!/usr/bin/env python
# encoding: utf-8
"""
plate_utils.py

Created by José Sánchez-Gallego on 28 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    28 Aug 2014 J. Sánchez-Gallego
      Initial version
   3 May 2015 J. Sánchez-Gallego
      Major rewrite

"""

from __future__ import division, print_function

import collections
import itertools

import numpy as np
from scipy.misc import factorial

from Totoro import config, exceptions, log, site
from Totoro.db import getConnection
from Totoro.utils import intervals


def updatePlate(plate, rearrangeIncomplete=False, **kwargs):
    """Finds new exposures and assigns them a new set.

    If `rearrangeIncomplete=True`, exposures in incomplete sets are then
    arranged in the best possible mode.

    """

    # updatePlate will likely fail if the code is being run within an open,
    # external session. So, we check to make sure that's not the case
    #  checkOpenSession()

    unassignedExposures = getUnassignedExposures(plate)

    newExposures = [exp for exp in unassignedExposures if exp.isValid(force=True, flag=True)[0]]

    if len(newExposures) == 0:
        return False

    log.debug('plate_id={0}: found {1} new exposures'.format(plate.plate_id, len(newExposures)))

    for exp in newExposures:
        assignExposureToOptimalSet(plate, exp)

        if rearrangeIncomplete:
            result = rearrangeSets(plate, mode='optimal', scope='incomplete', silent=True)

            if not result:
                return result

    return True


def getUnassignedExposures(plate):
    """Returns exposures in `plate` that are not assigned to a set."""

    from Totoro.dbclasses import Exposure as TotoroExposure

    scienceExposures = plate.getScienceExposures()

    unassigned = [
        TotoroExposure(exp) for exp in scienceExposures if exp.mangadbExposure[0].set_pk is None
    ]

    unassignedSorted = sorted(unassigned, key=lambda exp: exp.exposure_no)

    return unassignedSorted


def assignExposureToOptimalSet(plate, exposure):
    """Assigns `exposure` to the best possible set in `plate`.

    Loops over all the sets in `plate` to find the best placement for a new
    `exposure`. If no set is suitable, creates a new one.

    """

    from Totoro.dbclasses import Set as TotoroSet

    db = plate.db
    session = plate.session

    optimalSet = getOptimalSet(plate, exposure)

    if optimalSet is None:
        setPK = int(getConsecutiveSets(1)[0])
        with session.begin():
            if session.query(db.mangaDB.Set).get(setPK) is None:
                newSet = db.mangaDB.Set(pk=setPK)
                session.add(newSet)
                session.flush()
                assert newSet.pk is not None, \
                    'something failed while creating a new set'
                exposure.mangadbExposure[0].set_pk = newSet.pk

                log.debug('plate_id={0}: exposure_no={1} assigned to '
                          'new set pk={2}'.format(plate.plate_id, exposure.exposure_no, newSet.pk))
                totoroNewSet = TotoroSet(newSet)
                plate.sets.append(totoroNewSet)
            else:
                log.debug('plate_id={0}: something failed while assigning new '
                          'set_pk to exposure_no={1}'.format(plate.plate_id, exposure.exposure_no))
                return
    else:
        with session.begin():
            exposure.mangadbExposure[0].set_pk = optimalSet.pk
        for ss in plate.sets:
            if ss.pk == optimalSet.pk:
                ss.totoroExposures.append(exposure)

        log.debug('plate_id={0}: exposure_no={1} assigned to set pk={2}'.format(
            plate.plate_id, exposure.exposure_no, optimalSet.pk))

    return


def getOptimalSet(plate, exposure):
    """Returns the best set in `plate` for an `exposure` or None."""

    from Totoro.dbclasses import Set as TotoroSet

    dither = exposure.ditherPosition

    incompleteSets = [
        set for set in plate.sets if set.getStatus()[0] in ['Incomplete', 'Unplugged']
    ]

    validSets = []
    signalNoise = []
    for ss in incompleteSets:

        setDithers = ss.getDitherPositions()

        if dither in setDithers:
            # If the dither exists, skips this set
            continue
        elif dither is None:
            # If the exposure has not a dither position (usually for mock
            # mock exposures), selects one of the unused dithers in the set.
            tmpDither = getValidDither(ss)
            assert tmpDither is not None, 'failed getting valid dither'
            exposure.ditherPosition = tmpDither

        exposures = ss.totoroExposures + [exposure]
        mockSet = TotoroSet.fromExposures(exposures)
        status = mockSet.getStatus(silent=True)[0]

        if status in ['Good', 'Excellent']:
            validSets.append(ss)
            # Adds 100 to SN2 array to make sure complete set are always chosen
            signalNoise.append(mockSet.getSN2Array() + 100)
        elif status in ['Incomplete', 'Unplugged']:
            validSets.append(ss)
            signalNoise.append(mockSet.getSN2Array())

    # Restore original dither position, in case we have changed it
    exposure.ditherPosition = dither

    if len(validSets) == 0:
        return None

    signalNoise = np.array(signalNoise)

    # Calculates the contribution of each mock set to the total completion.
    completion = np.zeros((signalNoise.shape[0], 2), np.float)
    completion[:, 0] = np.nanmean(signalNoise[:, 0:2], axis=1)
    completion[:, 0] /= config['SN2thresholds']['plateBlue']
    completion[:, 1] = np.nanmean(signalNoise[:, 2:], axis=1)
    completion[:, 1] /= config['SN2thresholds']['plateRed']
    completion = np.nanmin(completion, axis=1)

    # Selects the set that contributes more to the total completion.
    return validSets[np.argmax(completion)]


def getValidDither(ss):
    """Returns a valid dither in a set."""

    ditherPositions = set(config['set']['ditherPositions'])
    setDithers = set(ss.getDitherPositions())

    if len(setDithers) == len(ditherPositions):
        return None

    remainingDithers = list(ditherPositions - setDithers)

    return remainingDithers[0]


def _getSetStatusLabel(exposure):
    """Returns the set status for an exposure or None."""

    try:
        return exposure._mangaExposure.set.status.label
    except (AttributeError, IndexError):
        return None


def rearrangeSets(plate, mode='complete', scope='all', force=False, LST=None, silent=False):
    """Rearranges exposures in a plate.

    If `mode='complete'`, uses a brute-force approach to obtain the best
    possible arrangement for exposures into sets. If `mode='sequential'`,
    removes the set allocation for the exposures in the plate and adds them
    back in the order in which they were observed.

    In `mode='complete'`, the function determines all the possible combinations
    of exposures in sets and chooses the one that maximises the SN2. Of all the
    permutations with SN2 within `config.setRearrangement.factor * maxSN2` and
    `maxSN2`, selects the one that leaves most of the incomplete sets at the
    beginning of the visibility window of the plate or to `LST`.

    Sets overridden good or bad are not affected by the rearrangement.

    Parameters
    ----------
    plate : `Totoro.Plate` instance
        The plate to rearrange.
    mode : str
        The mode of the rearrangement, either `'complete'` or `'sequential'`.
    scope : str
        If `'all'`, all the exposures in the plate are rearranged. If
        `'incomplete'`, only exposures in incomplete sets are used.
    force : bool
        The limit of permutations to tests depends on the number of exposures,
        and can be calculated using `getNumberPermutations`.
        config.setArrangement.permutationLimitPlate and
        config.setArrangement.permutationLimitIncomplete determine the maximum
        number of permutations to test. If `force` is True, those limits are
        ignored.
    LST : None or float
        If multiple permutations have optimal SN2, the one that moves the
        average of the medium points of the incomplete sets closer to this LST
        is selected. If `LST=None`, the current LST is used.
    silent : bool
        Determines the verbosity of the function

    Returns
    -------
    result : bool
        Returns True if the rearrangement was successful, False otherwise.

    """

    from Totoro.dbclasses import Exposure
    from Totoro.dbclasses import Set

    # Sets logging level
    if silent:
        logMode = log.debug
    else:
        logMode = log.info

    # Selects exposures to consider
    if scope.lower() == 'all':
        permutationLimit = config['setArrangement']['permutationLimitPlate']
        exposures = [Exposure(exp) for exp in plate.getScienceExposures()]
    elif scope.lower() == 'incomplete':
        permutationLimit = config['setArrangement']['permutationLimitIncomplete']
        exposures = []
        for ss in plate.sets:
            if ss.getStatus()[0] in ['Incomplete', 'Unplugged']:
                for exposure in ss.totoroExposures:
                    if not isinstance(exposure, Exposure):
                        exposure = Exposure(exposure)
                    exposures.append(exposure)
    else:
        raise exceptions.TotoroError('scope={0} is invalid'.format(scope))

    # Removes exposures that are in sets overriden good or bad, or that are
    # invalid.
    validExposures = []
    for exp in exposures:
        setStatus = _getSetStatusLabel(exp)
        if setStatus is not None and 'Override' in setStatus:
            continue
        elif not exp.isValid(force=True, flag=True)[0]:
            continue
        elif exp.isMock:
            validExposures.append(exp)
        else:
            validExposures.append(exp)

    # Stores overridden sets
    overridenSets = [
        ss for ss in plate.sets if ss.status is not None and 'Override' in ss.status.label
    ]

    # Does some logging.
    logMode('plate_id={0}: rearranging sets, mode=\'{1}\', scope=\'{2}\''.format(
        plate.plate_id, mode, scope))
    logMode('plate_id={0}: found {1} valid exposures'.format(plate.plate_id, len(validExposures)))

    # No exposures to consider.
    if len(validExposures) == 0:
        return True

    session = plate.session

    if mode.lower() == 'sequential':
        # If mode is sequential, removes set_pk from all selected exposures
        # and triggers a plate update.
        with session.begin():
            for exposure in validExposures:
                if exposure.mangadbExposure[0].set_pk is not None:
                    session.delete(exposure.mangadbExposure[0].set)
                exposure.mangadbExposure[0].set_pk = None
                exposure.mangadbExposure[0].exposure_status_pk = None

        plate.sets = []
        updatePlate(plate, rearrangeIncomplete=False)

        return True

    # The remainder of the function assumes that the mode is optimal.

    ditherPositions = [exp.ditherPosition for exp in validExposures]
    nPermutations = getNumberPermutations(ditherPositions)

    logMode('plate_id={0}: testing {1} permutations'.format(plate.plate_id, nPermutations))

    if nPermutations > permutationLimit:
        if force is False:
            logMode('plate_id={0}: hard limit for number of permutations '
                    'in rearrangement ({1}) reached. Not rearranging.'.format(
                        plate.plate_id, permutationLimit))
            return False
        else:
            logMode('plate_id={0}: hard limit for number of permutations '
                    'in rearrangement reached but ignoring because '
                    'force=True.'.format(plate.plate_id))

    permutations = calculatePermutations(ditherPositions)

    def getSetId(ss):
        """Returns unique identifier for a set from the ids of its exposures"""
        return np.sum([id(exp) for exp in ss.totoroExposures])

    zeroSN2 = np.array([0.0, 0.0, 0.0, 0.0])

    goodArrangements = []
    setStatus = {}
    setSN2 = {}
    completions = []

    # Adds the SN2 of the overridden sets
    for ss in overridenSets:
        setId = getSetId(ss)
        setStatus[setId] = ss.getStatus(silent=True)[0]
        if 'Good' in setStatus[setId]:
            setSN2[setId] = ss.getSN2Array()
        else:
            setSN2[setId] = zeroSN2

    setRearrFactor = config['setArrangement']['factor']

    for nn, permutation in enumerate(permutations):

        sets = []

        for setIndices in permutation:

            setExposures = [validExposures[ii] for ii in setIndices if ii is not None]

            ss = Set.fromExposures(setExposures)
            sets.append(ss)

            # To avoid calculating the state of a set more than one, creates
            # a dictionary with the quality of the set
            setId = getSetId(ss)
            if setId not in setStatus:
                setStatus[setId] = ss.getStatus(silent=True)[0]
                setSN2[setId] = ss.getSN2Array() \
                    if setStatus[setId] in ['Excellent', 'Good',
                                            'Override Good'] else zeroSN2

            del ss

        sets += overridenSets

        # Instead of using Plate.getPlateCompletion, we calculate the plate
        # completion here using the setStatus dictionary. Way faster this way.
        plateSN2 = np.nansum([setSN2[getSetId(ss)] for ss in sets], axis=0)

        blueSN2 = np.nanmean(plateSN2[0:2])
        blueCompletion = blueSN2 / config['SN2thresholds']['plateBlue']
        redSN2 = np.nanmean(plateSN2[2:])
        redCompletion = redSN2 / config['SN2thresholds']['plateRed']
        plateCompletion = np.min([blueCompletion, redCompletion])

        # If the plate completion is lower than setRearrFactor times the
        # current maximum completion, we don't bother storing this permutation.
        if (len(completions) == 0 or plateCompletion >= setRearrFactor * np.max(completions)):
            completions.append(plateCompletion)
            goodArrangements.append(fixBadSets(sets))

        # Every 10% of the permutations.
        if (nn + 1) * 100. / nPermutations % 10 == 0:
            logMode('{0:d}% completed'.format(int((nn + 1) * 100. / nPermutations)))

    logMode('{0} permutations tested.'.format(nPermutations))

    completions = np.array(completions)

    # If the scope is 'incomplete', adds the completion of the good sets.
    if scope.lower() == 'incomplete':
        plateCompletion = plate.getPlateCompletion(includeIncompleteSets=False)
        completions += plateCompletion

    # From the good arrangements already selected, find the optimal one.
    optimalArrangement = selectOptimalArrangement(goodArrangements, completions, LST=LST)

    # If the scope is 'incomplete', adds the good sets to the optimal
    # arrangement.
    if scope == 'incomplete':
        optimalArrangement = list(optimalArrangement)
        for ss in plate.sets:
            if (ss.getStatus(silent=True)[0] in ['Good', 'Excellent'] and
                    ss not in optimalArrangement):
                optimalArrangement.append(ss)

    # Applies the new arrangement and modifies the plate info accordingly.
    status = applyArrangement(plate, optimalArrangement)

    return status


def selectOptimalArrangement(arrangements, completions, LST=None):
    """Selects the best possible option from a list of set arrangements."""

    arrangements = np.array(arrangements)

    # If one of the arrangements completes the plate, selects the one with the
    # largest completion and the fewest sets.
    if np.any(completions > 1):
        complete = arrangements[completions == np.max(completions)]
        if len(complete) == 1:
            return complete[0]
        else:
            nSets = np.array([len(sets) for sets in complete])
            completeMinSets = complete[np.argmin(nSets)]
            return completeMinSets

    # If no complete plates exist, divides the completion by the number of sets
    nSets = np.array([len(sets) for sets in arrangements])
    completions = completions / nSets

    # Selects the top tier arrangements.
    setRearrFactor = config['setArrangement']['factor']
    minCompletion = np.max(completions) * setRearrFactor

    topArrangements = arrangements[completions >= minCompletion]

    # If only one arrangement exists, we are done
    if len(topArrangements) == 1:
        return topArrangements[0]

    # If several top arrangements exist, we select the one that has more
    # incomplete sets after the selected LST.
    if LST is None:
        LST = site.localSiderealTime()

    # For each arrangements, calculates the difference between `LST` and the
    # middle point of the LST window for each sets and sums them.
    cumulatedLSTdiffs = []
    for arrangement in topArrangements:
        LSTdiff = []
        for ss in arrangement:
            setMeanLST = intervals.calculateMean(ss.getLST(), wrapAt=24.)
            LSTdiff.append((setMeanLST - LST) % 24.)
        cumulatedLSTdiffs.append(np.sum(LSTdiff))

    cumulatedLSTdiffs = np.array(cumulatedLSTdiffs)

    # Returns the arrangement with the smallest cumulated LST
    return topArrangements[np.argmin(cumulatedLSTdiffs)]


def applyArrangement(plate, arrangement):
    """Applies a set arrangement to a plate."""

    from Totoro.dbclasses import Set

    db = plate.db
    session = plate.session

    # If all exposures are real, saves data to the DB.
    expMock = [exp.isMock for ss in arrangement for exp in ss.totoroExposures]

    if not any(expMock):

        # Selects only new sets and skips overridden sets
        arrangement = [ss for ss in arrangement
                       if ss.status is None or 'Override' not in ss.status.label]

        # Removes sets and exposure-set assignment from the DB
        with session.begin():
            for ss in plate.sets:
                if ss.status is not None and 'Override' in ss.status.label:
                    continue
                for exp in ss.totoroExposures:
                    setPK = exp.mangadbExposure[0].set_pk
                    exp.mangadbExposure.set_pk = None
                    if setPK is not None:
                        setDB = session.query(db.mangaDB.Set).get(setPK)
                        if setDB is not None:
                            session.delete(setDB)
                            session.flush()

                session.flush()

        # Gets the pks to use
        pks = getConsecutiveSets(len(arrangement))

        # Now creates the new sets and assigns the exposures
        with session.begin():
            for ii, ss in enumerate(arrangement):
                newSet = db.mangaDB.Set(pk=int(pks[ii]))
                session.add(newSet)
                session.flush()
                for exp in ss.totoroExposures:
                    exp.mangadbExposure[0].set_pk = newSet.pk

                    log.debug('plate_id={0}: exposure_no={1} assigned '
                              'to set pk={2}'.format(plate.plate_id, exp.exposure_no, newSet.pk))

        # Finally, reloads the exposures and sets into plate.
        plate.sets = []
        for ss in plate.getMangaDBSets():
            plate.sets.append(Set(ss))

    else:

        # If any of the exposures or sets is mock, just updates the data in the
        # plate, but saves nothing to the DB.

        plate.sets = arrangement

    return True


def calculatePermutations(inputList):
    """Calculates all the possible permutations of a list of dithers."""

    pairs = [(nn, inputList[nn]) for nn in range(len(inputList))]
    pairs = sorted(pairs, key=lambda value: value[1])

    splitPairs = [list(bb) for aa, bb in itertools.groupby(pairs, lambda value: value[1])]
    sortedPairs = sorted(splitPairs, key=lambda xx: len(xx))[::-1]

    indicesSeed = [[element[0] for element in sP] for sP in sortedPairs]
    nExpPerDither = [len(ii) for ii in indicesSeed]
    for ii in indicesSeed:
        if len(ii) < np.max(nExpPerDither):
            ii += [None] * (np.max(nExpPerDither) - len(ii))

    if len(indicesSeed) > 0:
        indices = [[tuple(indicesSeed[0])]]
        indices += [list(itertools.permutations(idx)) for idx in indicesSeed[1:]]
    else:
        indices = []

    cartesianProduct = itertools.product(*indices)
    for prod in cartesianProduct:
        yield list(itertools.zip_longest(*prod))


def getNumberPermutations(ditherPositions):
    """Estimates the number of permutations for a list of dithers."""

    if len(ditherPositions) == 0:
        return 0

    repDict = collections.defaultdict(int)
    for cc in ditherPositions:
        repDict[cc] += 1

    maxNDither = 0
    for key in list(repDict.keys()):
        if repDict[key] > maxNDither:
            maxNDither = repDict[key]

    return int(factorial(maxNDither)**(len(list(repDict.keys())) - 1))


def fixBadSets(sets):
    """Splits bad sets into a series of valid, incomplete sets."""

    from Totoro.dbclasses import Set as Set

    toRemove = []
    toAdd = []

    for ss in sets:

        if ss.getStatus(silent=True)[0] != 'Bad':
            continue

        toRemove.append(ss)

        if len(ss.totoroExposures) == 1:
            raise exceptions.TotoroError(
                'found bad set with one exposure. This is probably a bug.')
        elif len(ss.totoroExposures) == 2:
            # If the bad set has two exposures, splits it.
            toAdd += [Set.fromExposures(exp) for exp in ss.totoroExposures]
        else:
            # Tests all possible combinations of two exposures to check if one
            # of them is a valid set.
            validSets = []
            for ii, jj in [[0, 1], [0, 2], [1, 2]]:
                testSet = Set.fromExposures([ss.totoroExposures[ii], ss.totoroExposures[jj]])
                if testSet.getStatus(silent=True)[0] != 'Bad':
                    validSets.append(testSet)

            if len(validSets) == 0:
                # If no valid combinations, each exposures goes to a set.
                toAdd += [Set.fromExposures(exp) for exp in ss.totoroExposures]
            else:
                # Otherwise, selects the combination that produces an
                # incomplete set with maximum SN2.
                signalToNoise = [np.nansum(xx.getSN2Array()) for xx in validSets]

                maxSet = validSets[np.argmax(signalToNoise)]

                toAdd.append(maxSet)
                missingExposure = [
                    exp for exp in ss.totoroExposures if exp not in maxSet.totoroExposures
                ]

                toAdd.append(Set.fromExposures(missingExposure))

    for ss in toRemove:
        sets.remove(ss)

    for ss in toAdd:
        sets.append(ss)

    return sets


def getConsecutiveSets(nSets=1):
    """Returns a list of consecutive set pks that are not assigned."""

    db = getConnection()
    session = db.Session()

    # Finds already used set pks
    with session.begin():
        setPKs = session.query(db.mangaDB.Set.pk).all()

    # Creates a list of unused set pks
    setPKs = np.array(setPKs).squeeze().tolist()
    setPKs = sorted(setPKs)
    candidatePKs = np.array([ii for ii in range(1, setPKs[-1] + 1) if ii not in setPKs])

    # Splits the pks into groups of consecutive values
    candidatePKsSplit = np.split(candidatePKs, np.where(np.diff(candidatePKs) != 1)[0] + 1)

    # If there is a groups with at least as many values as nSets, uses it.
    pks = None
    for split in candidatePKsSplit:
        if len(split) >= nSets:
            pks = split[0:nSets]
            break

    # If no consecutive range of pks exists, just continues from the last pk
    if pks is None:
        pks = list(range(setPKs[-1] + 1, setPKs[-1] + 1 + nSets))

    return pks


def removeOrphanedSets():
    """Removes sets without exposures."""

    nRemoved = 0

    db = getConnection()
    session = db.Session()

    with session.begin():
        sets = session.query(db.mangaDB.Set).all()
        for ss in sets:
            if len(ss.exposures) == 0:
                session.delete(ss)
                nRemoved += 1

    log.debug('removed {0} orphaned sets'.format(nRemoved))

    return
