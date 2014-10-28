#!/usr/bin/env python
# encoding: utf-8
"""
mangaLogicTimeline.py

Created by José Sánchez-Gallego on 3 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    3 Aug 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import numpy as np
from sdss.internal.manga.Totoro import utils
from sdss.internal.manga.Totoro import config, site


expTime = (config['exposure']['exposureTime'] /
           config['plugger']['efficiency'])
maxAlt = config['exposure']['maxAltitude']


def getOptimalPlate(plates, jdRanges, prioritisePlugged=True, **kwargs):
    """Gets the optimal plate to observe in a range of JDs."""

    jdRanges = np.atleast_2d(jdRanges).copy()
    incompletePlates = [plate for plate in plates if plate.isComplete is False]

    if prioritisePlugged:
        pluggedPlates = [plate for plate in incompletePlates
                         if plate.isPlugged]
        if len(pluggedPlates) > 0:
            observedFlag = simulatePlates(pluggedPlates, jdRanges)
            if observedFlag is True:
                optimal = selectOptimal(pluggedPlates)
                cleanupPlates(pluggedPlates, optimal)
                return optimal

    notPlugged = [plate for plate in incompletePlates if not plate.isPlugged]

    startedPlates = [plate for plate in notPlugged
                     if plate.getPlateCompletion() > 0]
    if len(startedPlates) > 0:
        observedFlag = simulatePlates(startedPlates, jdRanges)
        if observedFlag is True:
            optimal = selectOptimal(startedPlates)
            cleanupPlates(startedPlates, optimal)
            return optimal

    observedFlag = simulatePlates(notPlugged, jdRanges)
    if observedFlag is False:
        return None

    optimal = selectOptimal(notPlugged)
    cleanupPlates(notPlugged, optimal)

    return optimal


def selectOptimal(plates, **kwargs):
    """Returns the optimal plate to observe."""

    plates = np.array(plates)

    newPlates = []
    for plate in plates:
        exps = plate.getTotoroExposures()
        for exp in exps:
            if hasattr(exp, '_tmp') and exp._tmp is True:
                newPlates.append(plate)
                break
    newPlates = np.array(newPlates)

    # First tries to select the complete plate with the fewer number
    # of exposures
    completedPlates = [plate for plate in newPlates if plate.isComplete]
    if len(completedPlates) > 0:
        nExposures = [len(plate.getTotoroExposures()) / plate.priority
                      for plate in completedPlates]
        # print('completed')
        return completedPlates[np.argmin(nExposures)]

    # If no completed plates, calculates the plate completion using only
    # complete sets
    plateCompletion = np.array(
        [plate.getPlateCompletion() for plate in newPlates])
    maxPlateCompletion = plateCompletion.max()

    # Selects the plates with maximum completion
    platesMaxCompletion = newPlates[np.where(
        plateCompletion == maxPlateCompletion)]

    if len(platesMaxCompletion) == 1:
        # If only one plate with maximum completion
        return platesMaxCompletion[0]
    else:
        # If several plates with maximum completion, returns that with the
        # maximum overall completion (including incomplete sets).
        plateCompletionIncomplete = np.array(
            [plate.getPlateCompletion(includeIncompleteSets=True)
             for plate in platesMaxCompletion])
        return platesMaxCompletion[plateCompletionIncomplete.argmax()]

    return None


def simulatePlates(plates, jdRanges):
    """Simulates exposures for a list of plates withing a range of JDs."""

    jdRanges = np.atleast_2d(jdRanges)
    observedFlag = False

    for plate in plates:
        plateLST = plate.getLSTRange()
        for jdRange in jdRanges:
            if not isObservable(plate, jdRange):
                continue
            jd = jdRange[0]
            while jd < jdRange[1]:

                if plate.isComplete and observedFlag is True:
                    return observedFlag

                lst = site.localSiderealTime(jd)
                if (not utils.isPointInInterval(lst, plateLST, wrapAt=24)):
                    pass
                elif plate.getAltitude(lst) > maxAlt:
                    pass
                else:
                    result = plate.addMockExposure(
                        set=None, startTime=jd, expTime=expTime, silent=True)
                    if result is not False:
                        observedFlag = True
                        result._tmp = True

                jd += expTime / 86400.

    return observedFlag  # True if we have added at least one exposure


def isObservable(plate, jdRanges):
    """Returns True if a plate is observable in a range of JDs."""

    jdRanges = np.atleast_2d(jdRanges)

    for jdRange in jdRanges:
        plateLSTRange = plate.getLSTRange()
        lstRange = site.localSiderealTime(jdRange)

        if (utils.getIntervalIntersection(plateLSTRange, lstRange, wrapAt=24.)
                is not False):
            return True

    return False


def cleanupPlates(plates, optimumPlate):
    """Clears exposures in plates that are not the optimum and removes the
    temporary flag in the optimum one."""

    for plate in plates:

        if plate is optimumPlate:
            for exp in plate.getValidExposures():
                if hasattr(exp, '_tmp'):
                    delattr(exp, '_tmp')
            for set in plate.sets:
                if len(set.totoroExposures) == 0:
                    plate.sets.remove(set)
            continue

        for set in plate.sets:
            set.totoroExposures = [exp for exp in set.totoroExposures
                                   if not hasattr(exp, '_tmp') or
                                   exp._tmp is not True]
        plate.sets = [set for set in plate.sets
                      if len(set.totoroExposures) > 0]

    return
