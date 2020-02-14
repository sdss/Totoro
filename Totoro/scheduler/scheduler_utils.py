#!/usr/bin/env python
# encoding: utf-8
"""
scheduler_utils.py

Created by José Sánchez-Gallego on 3 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    3 Aug 2014 J. Sánchez-Gallego
      Initial version
    6 May 2015 J. Sánchez-Gallego
      Major changes to the logic and structure

"""

from __future__ import division, print_function

from collections import OrderedDict
from numbers import Number

import numpy as np

from Totoro import config, site, utils
from Totoro.scheduler import observingPlan


expTime = config['exposure']['exposureTime']
minSchedulingTime = config['plugger']['minSchedulingTimeHours']
patchSetFactor = config['scheduling']['patchSetFactor']
platePriorityFactor = config['scheduling']['platePriorityFactor']
nextNightFactor = config['scheduling']['nextNightFactor']


def _getNextNightRange(jdRange):
    """For a given jdRange, returns the JD range of the next night."""

    currentNight = observingPlan[(observingPlan['JD0'] <= jdRange[0]) &
                                 (observingPlan['JD1'] >= jdRange[1])]

    if len(currentNight) == 0 or len(currentNight) > 1:
        return None

    JD = currentNight['JD']

    # Next night must be JD + 1 (i.e., be in the same run)
    nextNight = observingPlan[observingPlan['JD'] == JD + 1]
    if len(nextNight) == 0:
        return None
    else:
        return (nextNight['JD0'][0], nextNight['JD1'][0])


def _addBookkeepingAttrs(plates):
    """Stores some information about the plates before the simulation"""

    for plate in plates:
        plate._before = {}
        plate._before['completion'] = plate.getPlateCompletion(useMock=True)
        plate._before['completion+'] = plate.getPlateCompletion(
            useMock=True, includeIncompleteSets=True)
        plate._before['nExposures'] = len(plate.getTotoroExposures(onlySets=True))


def getDictOfSchedulablePlates(plates, mode):
    """Returns an OrderedDict with categories of scheduleable plates.

    This function creates an ordered dictionary in which the keys are the
    categories of plates in the order in which they should be scheduled to
    look for an optimal plate.

    If `mode='plugger'` the categories are `'plugged'`, `'platesWithSignal'`,
    `'incomplete'`, and `'backup'`.

    For `mode='planner'`, categories are `'platesWithSignal'`,
    `'drilled'`, `'fieldsInFootprint'`, `'backup'`, and
    `'fieldsOutsideFootprint'`.

    """

    from Totoro.dbclasses.plate import Plate
    from Totoro.dbclasses.field import Field

    assert mode in ['planner', 'plugger']

    if mode == 'plugger':
        schPlates = OrderedDict([('plugged', []), ('platesWithSignal', []), ('incomplete', []),
                                 ('backup', [])])
        for plate in plates:
            if plate.isPlugged:
                schPlates['plugged'].append(plate)
            elif plate.getPlateCompletion(includeIncompleteSets=True) > 0:
                schPlates['platesWithSignal'].append(plate)
            elif plate.statuses[0].label.lower() != 'backup':
                schPlates['incomplete'].append(plate)
            else:
                schPlates['backup'].append(plate)

    elif mode == 'planner':
        schPlates = OrderedDict([('platesWithSignal', []),
                                 ('drilled', []),
                                 ('deep_plates', []),
                                 ('fieldsInFootprint', []),
                                 ('fieldsOutsideFootprint', []),
                                 ('backup', []),
                                 ('backup_outside', [])])
        for plate in plates:
            isPlate = isinstance(plate, Plate) and not isinstance(plate, Field)
            isBackup = (hasattr(plate, 'statuses') and len(plate.statuses) > 0 and
                        plate.statuses[0].label.lower() == 'backup')
            isInFootPrint = plate.inFootprint
            if (plate.getPlateCompletion(includeIncompleteSets=True) > 0 and
                    plate.completion_factor <= 1):
                schPlates['platesWithSignal'].append(plate)
            elif isPlate and plate.completion_factor <= 1:
                schPlates['drilled'].append(plate)
            elif isPlate and plate.completion_factor > 1:
                schPlates['deep_plates'].append(plate)
            # elif ((isPlate and not isBackup) or
            #       (plate.dateAtAPO is not None and plate.dateAtAPO > 0)):
            #     if isInFootPrint:
            #         schPlates['backup'].append(plate)
            #     else:
            #         schPlates['backup_outside'].append(plate)
            elif not isPlate and isInFootPrint:
                schPlates['fieldsInFootprint'].append(plate)
            elif isPlate and isBackup:
                if isInFootPrint:
                    schPlates['backup'].append(plate)
                else:
                    schPlates['backup_outside'].append(plate)
            elif not isPlate and not isInFootPrint:
                schPlates['fieldsOutsideFootprint'].append(plate)

    # Performs a couple sanity checks. Makes sure each plate is in one and
    # and only one category.
    allPlatesFromDict = [plate for cc in schPlates for plate in schPlates[cc]]
    assert len(allPlatesFromDict) == len(set(allPlatesFromDict)), \
        'there are plates in more than one scheduling category.'
    # assert set(allPlatesFromDict) == set(plates), \
    #     'some plates are not being assigned a scheduling category.'

    return schPlates


def getOptimalPlate(plates, jdRange, mode='plugger', prioritiseAPO=None,
                    observedPlates=[], **kwargs):
    """Gets the optimal plate to observe in a range of JDs."""

    assert len(jdRange) == 2

    if prioritiseAPO is None:
        prioritiseAPO = True if mode == 'plugger' else False

    optimalPlate = None

    # Makes sure we are not dealing with any completed plates.
    incompletePlates = [
        plate for plate in plates
        if not utils.isPlateComplete(plate, write_apocomplete=False, mark_complete=False)
    ]

    # Selects plates that intersect with the observing window for at least
    # one exposure.
    observablePlates = [plate for plate in incompletePlates if isObservable(plate, jdRange)]

    # If there are no plates that meet those requirements, uses all the
    # incomplete plates
    if len(observablePlates) == 0:
        # log.warning('no plates found that intersect with the beginning of '
        #             'the observing window. Scheduling all available plates.',
        #             TotoroUserWarning)
        observablePlates = incompletePlates

    # Creates the dictionary of scheduling categories
    schPlatesDict = {}
    schPlatesDict['nightObserved'] = observedPlates
    schPlatesDict.update(getDictOfSchedulablePlates(observablePlates, mode))

    # Loops over the list of categories, scheduling each list of plates
    # and trying to find an optimal plate. If it gets it, returns the
    # plate. Otherwise, jumps to the next category.
    for plateCat in schPlatesDict:
        platesToSchedule = schPlatesDict[plateCat]

        # Defines the normalise and scope modes to use depending on the
        # category and mode.
        if mode == 'plugger' and plateCat == 'plugged':
            normalise = False
            scope = 'plugged'
        elif plateCat == 'platesWithSignal':
            scope = 'plugged'
            normalise = True
        else:
            normalise = True
            scope = 'all'

        if not prioritiseAPO:

            optimalPlate, newExposures = runSimulation(
                platesToSchedule, jdRange, mode=mode, scope=scope, normalise=normalise)

            if optimalPlate is not None:
                return optimalPlate, newExposures

        else:

            for locPlates in _getPlatesAtAPO(platesToSchedule):
                optimalPlate, newExposures = runSimulation(
                    locPlates, jdRange, mode=mode, scope=scope, normalise=normalise)

                if optimalPlate is not None:
                    return optimalPlate, newExposures

    # No optimal plate found.
    return None, []


def _getPlatesAtAPO(plates):
    """Returns lists of plates at APO and outside it."""

    atAPO = []
    notAtAPO = []

    for plate in plates:
        if plate.getLocation() == 'APO':
            atAPO.append(plate)
        else:
            notAtAPO.append(plate)

    return atAPO, notAtAPO


def runSimulation(plates, jdRange, mode='plugger', scope='all', normalise=True, **kwargs):
    """Runs the simulation for a subset of plates."""

    # Adss bookkeeping information
    _addBookkeepingAttrs(plates)

    simulatePlates(plates, jdRange, mode, **kwargs)
    optimalPlate = selectPlate(plates, jdRange, scope=scope, normalise=normalise)

    if optimalPlate:
        newExps = cleanupPlates(plates, optimalPlate, jdRange)
        return optimalPlate, newExps
    else:
        cleanupPlates(plates, None, jdRange)
        return None, []


def _normaliseWindowLength(plates, jdRange, factor=1.0, apply=True):
    """Calculates normalisation factors based on window lengths."""

    lstRange = site.localSiderealTime(jdRange)
    plateWindowLength = np.array([
        utils.getIntervalIntersectionLength(plate.getLSTRange(), lstRange, wrapAt=24.)
        for plate in plates
    ])

    if len(plateWindowLength[plateWindowLength > 0]) == 0:
        return

    # When called for the next night, some plates may have length 0. We reject
    # those for the calculation of the minimum length.
    minLength = np.min(plateWindowLength[plateWindowLength > 0])

    plateWindowLenghtNormalised = np.array((plateWindowLength / minLength))

    plateWindowLenghtNormalisedFactor = (factor * (plateWindowLenghtNormalised - 1.) + 1)

    # We restore the original factor for plates with window length 0.
    # To be improved at some point.
    plateWindowLenghtNormalisedFactor[np.where(plateWindowLength == 0)] = 1

    if apply:
        _completionFactor(plates, plateWindowLenghtNormalisedFactor)

    return plateWindowLenghtNormalisedFactor


def selectPlate(plates, jdRange, normalise=False, scope='all'):
    """From a list of simulated plates, returns the optimal one."""

    # Gets the JD range for the following night
    nextNightJDrange = _getNextNightRange(jdRange)

    # First we exclude plates without new exposures
    plates = [plate for plate in plates if plate._after['nNewExposures'] > 0]

    # Sorts plates by inverse plate completion.
    plates = sorted(plates, reverse=True, key=lambda plate: plate.getPlateCompletion()
                    if plate.getPlateCompletion() <= 1 else 1. / plate.getPlateCompletion())

    if len(plates) == 0:
        return None

    # If we are scheduling only plugged plates, we rather plug a new plate
    # unless we can observe a plugged plate at least for a whole set.
    availableTime = (jdRange[1] - jdRange[0]) * 24.
    completionIncrease = np.array(
        [plate._after['completion'] - plate._before['completion'] for plate in plates])

    # minSchedulingTime ensures that if the remaining time < length of a set,
    # we still use the plugged plates, if any.
    if scope == 'plugged':
        if (availableTime > minSchedulingTime and np.all(completionIncrease == 0)):
            return None
    else:
        # If no plate has been observed for a whole set, we try to use first
        # plates that are already plugged.
        if np.all(completionIncrease == 0):
            pluggedPlates = [plate for plate in plates if plate.isPlugged]
            if len(pluggedPlates) > 0:
                plates = pluggedPlates

    # If plugger, tries to select only plates at APO
    if scope == 'plugged':

        platesAtAPO = [plate for plate in plates if plate.getLocation() == 'APO']
        if len(platesAtAPO) > 0:
            plates = platesAtAPO

        # Now tries to select only plates that have been marked.
        markedPlates = [
            plate for plate in plates if 'Accepted' in [status.label for status in plate.statuses]
        ]
        if len(markedPlates) > 0:
            plates = markedPlates

    # We check if any of the plate is complete after the simulation.
    # If so, we return the one with fewer new exposures.
    completePlates = [plate for plate in plates
                      if plate._after['completion'] > plate.completion_factor]
    nNewExposures = [plate._after['nNewExposures'] for plate in completePlates]
    if len(completePlates) > 0:
        return completePlates[np.argmin(nNewExposures)]

    # We record the real completion before and after. We will normalise the
    # other completions based on our scheduling logic.
    for plate in plates:
        plate._before['realCompletion'] = plate._before['completion']
        plate._before['realCompletion+'] = plate._before['completion+']
        plate._after['realCompletion'] = plate._after['completion']
        plate._after['realCompletion+'] = plate._after['completion+']

    # If normalise=True, we divide the several completion values by the
    # length of the observing window for the plate, normalised by the length
    # of the minimum plate window. The effect of this is that plates with short
    # observing windows get comparatively larger completions and, thus, have
    # higher chance of being selected. This is good for plugged plates, as it
    # tries to schedule first plates with short windows even if other plates
    # could be completed at the time.

    # We also increase the completion of plates for which we have patched sets,
    # while we penalise those with incomplete sets. With this logic, we hope
    # that plates are observed when their incomplete sets can be patched.

    if normalise:

        _normaliseWindowLength(plates, jdRange, factor=1.0, apply=True)

        # We also normalise using the following night, if possible.
        if nextNightJDrange is not None:
            _normaliseWindowLength(plates, nextNightJDrange, factor=nextNightFactor, apply=True)

        # Now we normalise plate completion using a metric that gives higher
        # priority to plates for which we have patched incomplete sets.
        patchedSetFactor = []
        for plate in plates:
            nSetsFactor = 0
            for ss in plate.sets:
                if not ss.isMock:
                    nNewExps = 0
                    for exp in ss.totoroExposures:
                        if hasattr(exp, '_tmp') and exp._tmp:
                            nNewExps += 1
                    setComplete = ss.getStatus()[0] in ['Good', 'Excellent']
                    if setComplete and nNewExps == 0:
                        pass
                    else:
                        if nNewExps > 0:
                            nSetsFactor += 2 * nNewExps
                            if setComplete:
                                nSetsFactor *= 2
                        else:
                            nSetsFactor -= 1

            patchedSetFactor.append(1. + patchSetFactor * nSetsFactor)

        _completionFactor(plates, patchedSetFactor)

    # We add the priority into the mix
    platePriorities = np.array([plate.priority for plate in plates]) - 5.
    _completionFactor(plates, 1 + platePriorityFactor * platePriorities)

    ancillaryPriorities = []
    for plate in plates:
        if hasattr(plate, 'ancillary_weight'):
            ancillaryPriorities.append(plate.ancillary_weight)
        else:
            ancillaryPriorities.append(1)
    _completionFactor(plates, np.array(ancillaryPriorities))

    # Selects the plates that have the largest increase in completion
    completionIncrease = [plate._after['completion'] - plate._before['completion']
                          for plate in plates if plate.completion_factor <= 1.]

    if len(completionIncrease) == 0:
        for plate in plates:
            if plate.completion_factor > 1:
                completionIncrease.append(
                    plate._after['completion'] - plate._before['completion'])

    completionIncrease = np.array(completionIncrease)

    plates = np.array(plates)
    maxCompletionIncrease = np.max(completionIncrease)
    plates = plates[np.where(completionIncrease == maxCompletionIncrease)]

    if len(plates) == 1:
        return plates[0]

    # If maxCompletionIncrease is 0, it means that no plate has been
    # observed for at least a set. In this case, if possible, we want to use
    # a plate that already has signal.
    if maxCompletionIncrease == 0:
        platesWithSignal = [plate for plate in plates if plate._before['completion+'] > 0]
        if len(platesWithSignal) > 0:
            plates = platesWithSignal

    # If several plates have maximum completion increase, use the incomplete
    # sets to break the tie.
    completionIncreasePlus = np.array(
        [plate._after['completion+'] - plate._before['completion+'] for plate in plates])

    return plates[np.argmax(completionIncreasePlus)]


def _completionFactor(plates, factor):
    """Multiplies plate completion by a factor."""

    plates = np.atleast_1d(plates)
    if isinstance(factor, Number):
        factor = len(plates) * [factor]
    factor = np.array(factor)

    for ii, plate in enumerate(plates):
        plate._before['completion'] *= factor[ii]
        plate._before['completion+'] *= factor[ii]
        plate._after['completion'] *= factor[ii]
        plate._after['completion+'] *= factor[ii]

    return plates


def simulatePlates(plates,
                   jdRange,
                   mode,
                   efficiency=None,
                   SN2Factor=None,
                   maxAltitude=None,
                   **kwargs):
    """Simulates exposures for a list of plates within a range of JDs."""

    from .planner import IC342_range

    mode = mode.lower()
    efficiency = efficiency if efficiency else config[mode]['efficiency']
    SN2Factor = SN2Factor if SN2Factor else config[mode]['simulationFactor']
    maxAltitude = maxAltitude if maxAltitude else config[mode]['maxAltitude']

    # This flag let us know if at least one plate has been successfully
    # simulated.
    success = False

    # Defines if we are going to rearrange exposures in incomplete sets after
    # each new mock exposure is added to a plate
    rearrange = False  # True if mode == 'plugger' else False

    for plate in plates:

        plateLST = plate.getLSTRange()
        plugging = plate.getActivePlugging()

        initial_completion = plate.getPlateCompletion()

        jd = jdRange[0]
        while jd < jdRange[1]:

            expTimeEff = plate.exposure_time / efficiency

            # If MaNGA is observing at the beginning of the night the cart is
            # already loaded, so we can assume that the efficiency of the first
            # exposure is 100%.
            row = observingPlan[observingPlan['JD'] == int(jd)]
            if len(row) > 0:
                row = row[0]
                if row['Position'] == 1 and jd == row['JD0']:
                    expTimeEff = expTime

            if plate.getPlateCompletion() >= plate.completion_factor:
                break

            # For deep plates, only allow one complete set at a time
            if (plate.completion_factor > 1 and
                    plate.getPlateCompletion() > initial_completion):
                break

            # For IC342 check if plate is complete
            if plate.field_name == 'IC342' and utils.isPlateComplete(plate, mark_complete=False,
                                                                     write_apocomplete=False):
                break

            lst = site.localSiderealTime(jd)

            # LST reserved for IC342 plates
            if lst > IC342_range[0] and lst < IC342_range[1]:
                if plate.field_name != 'IC342':
                    break

            lstMean = site.localSiderealTime(jd + expTimeEff / 2. / 86400.)

            # Calculates how much time the exposure is observed above
            # max altitude minus 5 degrees.
            # lstRange = site.localSiderealTime([jd, jd + expTimeEff / 86400.])
            # highAltitudeRange = plate.getLSTRangeAboveAltitude(
            #     maxAltitude - 5)
            # highAltitudeIntersection = utils.getIntervalIntersectionLength(
            #     lstRange, highAltitudeRange, wrapAt=24.)

            if (not utils.isPointInInterval(lst, plateLST, wrapAt=24)):
                # if mode == 'planner':
                break
            elif plate.getAltitude(lstMean) > maxAltitude:
                break
            # elif highAltitudeIntersection * 3600 > expTimeEff * 0.9:
            #     # If the most part of the exposure happens at high altitude,
            #     # the exposure is not valid.
            #     break
            else:
                result = plate.addMockExposure(
                    set=None,
                    startTime=jd,
                    expTime=expTimeEff,
                    plugging=plugging,
                    factor=SN2Factor,
                    rearrange=rearrange,
                    silent=True)

                # If everything worked, we add a flag to mark this exposure
                # as temporary.
                if result is not False:
                    result._tmp = True
                    success = True
                else:
                    break

            jd += expTimeEff / 86400.

    # Now loops over the plates again and updates the dictionary with
    # simulation information.
    for plate in plates:
        plate._after = {}
        plate._after['completion'] = plate.getPlateCompletion(useMock=True)
        plate._after['completion+'] = plate.getPlateCompletion(
            useMock=True, includeIncompleteSets=True)
        plate._after['nNewExposures'] = len([
            exp for exp in plate.getTotoroExposures(onlySets=True)
            if hasattr(exp, '_tmp') and exp._tmp
        ])

    return success


def isObservable(plate, jdRange):
    """Returns True if a plate overlaps with the first element of jdRange
    and can be observed for at least one exposure."""

    plateLSTRange = plate.getLSTRange()
    lstRange = site.localSiderealTime(jdRange)

    if not utils.intervals.isPointInInterval(lstRange[0], plateLSTRange, wrapAt=24.):
        return False

    if utils.intervals.getIntervalIntersectionLength(
            plateLSTRange, lstRange, wrapAt=24.) < expTime / 3600.:
        return False

    return True


def cleanupPlates(plates, optimalPlate, jdRange):
    """Clears exposures for plates that are not optimal. Returns new exposures.

    For the optimal plate it removes exposures in incomplete sets
    if jRange > 1 hour.

    """

    # We start by removing everything from non-optimal plates
    for plate in plates:

        if plate is optimalPlate:
            continue

        for attr in ['_before', '_after', '_tmp']:
            if hasattr(plate, attr):
                delattr(plate, attr)

        for ss in plate.sets:
            ss.totoroExposures = [
                exp for exp in ss.totoroExposures if not hasattr(exp, '_tmp') or not exp._tmp
            ]
            ss._status = None

        plate.sets = [ss for ss in plate.sets if not ss.isMock or len(ss.totoroExposures) > 0]

    if optimalPlate is None:
        return

    # Calculates change in completion rate in the optimal plate
    completionChange = (optimalPlate._after['completion'] - optimalPlate._before['completion'])

    for attr in ['_before', '_after', '_tmp']:
        if hasattr(optimalPlate, attr):
            delattr(optimalPlate, attr)

    # Calculates how much time we have actually scheduled with the
    # optimal plate
    newExpJDs = np.array([
        exp.getJD() for exp in optimalPlate.getTotoroExposures()
        if hasattr(exp, '_tmp') and exp._tmp
    ])
    timeJDRange = np.sum(newExpJDs[:, 1] - newExpJDs[:, 0])

    # If the plate has a change in completion, we may want to exclude the new
    # exposures that do not form complete sets (that is, that do not contribute
    # to completion).
    if completionChange > 0:

        # Identifies the last set and check if it is incomplete.
        lastSet = optimalPlate.getLastSet()

        if lastSet.getStatus()[0] in ['Incomplete', 'Unplugged']:

            lastSet._status = None

            # Finds new exposures in the last set and calculates the jd range
            # they cover.
            exposuresToRemove = [
                exp for exp in lastSet.totoroExposures if hasattr(exp, '_tmp') and exp._tmp
            ]

            if len(exposuresToRemove) > 0:
                exposuresToRemoveJD = np.array([exp.getJD() for exp in exposuresToRemove])
                timeToRemove = np.sum(exposuresToRemoveJD[:, 1] - exposuresToRemoveJD[:, 0])

                # Calculates how much time unscheduled there is if we remove
                # the exposures.
                timeLeftAfterRemoval = (jdRange[1] - jdRange[0] - timeJDRange + timeToRemove)

                # If the time is > 1 hour, it might be that we can do something
                # more useful with that time, so we remove the exposures.
                if timeLeftAfterRemoval * 24. > minSchedulingTime:
                    for exp in exposuresToRemove:
                        lastSet.totoroExposures.remove(exp)

    # For the remaining mock exposures, we remove the _simulationData and _tmp
    # attributes.
    newExps = []
    for exp in optimalPlate.getTotoroExposures(onlySets=True):
        for attr in ['_before', '_after', '_tmp']:
            if hasattr(exp, attr):
                delattr(exp, attr)
                newExps.append(exp)

    # Makes sure there are no empty sets in the plate
    for ss in optimalPlate.sets:
        if ss.isMock and len(ss.totoroExposures) == 0:
            optimalPlate.sets.remove(ss)

    return newExps
