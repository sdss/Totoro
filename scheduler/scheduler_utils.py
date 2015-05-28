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

from __future__ import division
from __future__ import print_function
import numpy as np
from sdss.internal.manga.Totoro import utils
from sdss.internal.manga.Totoro import config, site
from sdss.internal.manga.Totoro.scheduler import observingPlan
from numbers import Number


expTime = config['exposure']['exposureTime']
minSchedulingTime = config['plugger']['minSchedulingTimeHours']


def getOptimalPlate(plates, jdRange, prioritisePlugged=True,
                    mode='plugger', **kwargs):
    """Gets the optimal plate to observe in a range of JDs."""

    assert len(jdRange) == 2

    optimalPlate = None

    # Makes sure we are not dealing with completed plates. Instead of using
    # plate.isComplete uses the plate completion so that it can account for
    # mock exposures (useful when the function is called with plates already)
    # simulated.
    incompletePlates = [plate for plate in plates
                        if plate.getPlateCompletion(useMock=True) < 1.]

    # Selects plates that overlap with the beginning of the JD range
    # and intersect with the observing window for at least one exposure.
    observablePlates = [plate for plate in incompletePlates
                        if isObservable(plate, jdRange)]

    # If there are no plates that meet those requirements, uses all the
    # incomplete plates
    if len(observablePlates) == 0:
        observablePlates = incompletePlates

    # Stores some information about the plates before the simulation
    for plate in observablePlates:
        plate._before = {}
        plate._before['completion'] = plate.getPlateCompletion(useMock=True)
        plate._before['completion+'] = plate.getPlateCompletion(
            useMock=True, includeIncompleteSets=True)
        plate._before['nExposures'] = len(
            plate.getTotoroExposures(onlySets=True))

    # Now selects the optimal plate
    # If mode is plugger. We try using the plates that are plugged. Note that
    # we use normalise=True here because we want to take into account the
    # observing windows of the plates already plugged.
    if mode == 'plugger':
        pluggedPlates = [plate for plate in observablePlates
                         if plate.isPlugged]

        if len(pluggedPlates) > 0:
            simulatePlates(pluggedPlates, jdRange, mode, **kwargs)
            optimalPlate = selectPlate(pluggedPlates, jdRange, scope='plugged',
                                       normalise=True)

            if optimalPlate:
                newExps = cleanupPlates(observablePlates, optimalPlate,
                                        jdRange)
                return optimalPlate, newExps
            else:
                cleanupPlates(pluggedPlates, None, jdRange)
                observablePlates = [plate for plate in observablePlates
                                    if plate not in pluggedPlates]

    # If that is not the case.
    simulatePlates(observablePlates, jdRange, mode, **kwargs)
    optimalPlate = selectPlate(observablePlates, jdRange)

    # Cleans up all plates
    newExps = cleanupPlates(observablePlates, optimalPlate, jdRange)

    return optimalPlate, newExps


def selectPlate(plates, jdRange, normalise=False, scope='all'):
    """From a list of simulated plates, returns the optimal one."""

    # First we exclude plates without new exposures
    plates = [plate for plate in plates
              if plate._after['nNewExposures'] > 0]

    if len(plates) == 0:
        return None

    # If we are scheduling only plugged plates, we rather plug a new plate
    # unless we can observe a plugged plate at least for a whole set.

    availableTime = (jdRange[1] - jdRange[0]) * 24.
    completionIncrease = np.array([plate._after['completion'] -
                                   plate._before['completion']
                                   for plate in plates])

    if scope == 'plugged' and availableTime > minSchedulingTime:
        if np.all(completionIncrease == 0):
            return None

    # Tries to select only plates at APO
    platesAtAPO = [plate for plate in plates if plate.getLocation() == 'APO']
    if len(platesAtAPO) > 0:
        plates = platesAtAPO

    # Now tries to select only plates that have been marked.
    markedPlates = [plate for plate in plates
                    if 'Accepted' in [status.label
                                      for status in plate.statuses]]
    if len(markedPlates) > 0:
        plates = markedPlates

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

        lstRange = site.localSiderealTime(jdRange)
        plateWindowLength = np.array([
            utils.getIntervalIntersectionLength(plate.getLSTRange(),
                                                lstRange, wrapAt=24.)
            for plate in plates])

        plateWindowLenghtNormalised = (plateWindowLength /
                                       np.min(plateWindowLength))

        _completionFactor(plates, 1.0 / plateWindowLenghtNormalised)

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

            patchedSetFactor.append(1. + 0.1 * nSetsFactor)

        _completionFactor(plates, patchedSetFactor)

    # We check if any of the plate is complete after the simulation.
    # If so, we return the one with fewer new exposures.
    completePlates = [plate for plate in plates
                      if plate._after['completion'] > 1]
    nNewExposures = [plate._after['nNewExposures'] for plate in completePlates]
    if len(completePlates) > 0:
        return completePlates[np.argmin(nNewExposures)]

    # We add the priority into the mix
    platePriorities = np.array([plate.priority for plate in plates])
    _completionFactor(plates, 1 + 0.25 * platePriorities)

    # If no complete plates exist, selects the ones that have the largest
    # increase in completion
    completionIncrease = np.array(
        [plate._after['completion'] - plate._before['completion']
         for plate in plates])

    plates = np.array(plates)

    platesMaxCompletionIncrease = plates[
        np.where(completionIncrease == np.max(completionIncrease))]

    if len(platesMaxCompletionIncrease) == 1:
        return platesMaxCompletionIncrease[0]

    # If several plates have maximum completion increase, use the incomplete
    # sets to break the tie
    completionIncreaseIncomplete = np.array(
        [plate._after['completion+'] - plate._before['completion+']
         for plate in platesMaxCompletionIncrease])

    optimalPlateWithIncomplete = platesMaxCompletionIncrease[
        np.argmax(completionIncreaseIncomplete)]

    return optimalPlateWithIncomplete


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


def simulatePlates(plates, jdRange, mode, efficiency=None, SN2Factor=None,
                   maxAltitude=None, **kwargs):
    """Simulates exposures for a list of plates within a range of JDs."""

    mode = mode.lower()
    efficiency = efficiency if efficiency else config[mode]['efficiency']
    SN2Factor = SN2Factor if SN2Factor else config[mode]['simulationFactor']
    maxAltitude = maxAltitude if maxAltitude else config[mode]['maxAltitude']

    # This flag let us know if at least one plate has been successfully
    # simulated.
    success = False

    for plate in plates:

        plateLST = plate.getLSTRange()
        plugging = plate.getActivePlugging()

        jd = jdRange[0]
        while jd < jdRange[1]:

            # If MaNGA is observing at the beginning of the night the cart is
            # already loaded, so we can assume that the efficiency of the first
            # exposure is 100%.
            row = observingPlan[observingPlan['JD'] == int(jd)]
            if len(row) > 0:
                row = row[0]
                if row['Position'] == 1 and jd == row['JD0']:
                    expTimeEff = expTime
                else:
                    expTimeEff = expTime / efficiency

            if plate.isComplete:
                break

            lst = site.localSiderealTime(jd)
            lstMean = site.localSiderealTime(jd+expTimeEff/2./86400.)

            if (not utils.isPointInInterval(lst, plateLST, wrapAt=24)):
                if mode == 'planner':
                    break
            elif plate.getAltitude(lstMean) > maxAltitude:
                break
            else:
                result = plate.addMockExposure(set=None, startTime=jd,
                                               expTime=expTimeEff,
                                               plugging=plugging,
                                               factor=SN2Factor, silent=True)

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
        plate._after['nNewExposures'] = len(
            [exp for exp in plate.getTotoroExposures(onlySets=True)
             if hasattr(exp, '_tmp') and exp._tmp])

    return success


def isObservable(plate, jdRange):
    """Returns True if a plate overlaps with the first element of jdRange
    and can be observed for at least one exposure."""

    plateLSTRange = plate.getLSTRange()
    lstRange = site.localSiderealTime(jdRange)

    if not utils.intervals.isPointInInterval(lstRange[0], plateLSTRange,
                                             wrapAt=24.):
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
        for ss in plate.sets:
            ss.totoroExposures = [exp for exp in ss.totoroExposures
                                  if not hasattr(exp, '_tmp') or not exp._tmp]
        plate.sets = [ss for ss in plate.sets
                      if not ss.isMock or len(ss.totoroExposures) > 0]

    if optimalPlate is None:
        return

    # Calculates change in completion rate
    completionChange = (optimalPlate._after['completion'] -
                        optimalPlate._before['completion'])

    # Calculates how much time we have actually scheduled with the
    # optimal plate
    newExpJDs = np.array(
        [exp.getJD() for exp in optimalPlate.getTotoroExposures()
         if hasattr(exp, '_tmp') and exp._tmp])
    timeJDRange = np.sum(newExpJDs[:, 1] - newExpJDs[:, 0])

    # If the plate has a change in completion, we may want to exclude the new
    # exposures that do not form complete sets (that is, that do not contribute
    # to completion).
    if completionChange > 0:

        # Identifies the last set and check if it is incomplete.
        lastSet = optimalPlate.getLastSet()

        if lastSet.getStatus()[0] in ['Incomplete', 'Unplugged']:

            # Finds new exposures in the last set and calculates the jd range
            # they cover.
            exposuresToRemove = [exp for exp in lastSet.totoroExposures
                                 if hasattr(exp, '_tmp') and exp._tmp]

            exposuresToRemoveJD = np.array(
                [exp.getJD() for exp in exposuresToRemove])
            timeToRemove = np.sum(exposuresToRemoveJD[:, 1] -
                                  exposuresToRemoveJD[:, 0])

            # Calculates how much time unscheduled there is if we remove the
            # exposures.
            timeLeftAfterRemoval = (jdRange[1] - jdRange[0] - timeJDRange +
                                    timeToRemove)

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
