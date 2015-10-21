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
from sdss.internal.manga.Totoro.exceptions import TotoroError


expTime = (config['exposure']['exposureTime'] /
           config['plugger']['efficiency'])
maxAlt = config['exposure']['maxAltitude']


def getOptimalPlate(plates, jdRanges, prioritisePlugged=True, **kwargs):
    """Gets the optimal plate to observe in a range of JDs."""

    jdRanges = np.atleast_2d(jdRanges)
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
    observedFlag = simulatePlates(notPlugged, jdRanges)
    if observedFlag is False:
        return None

    optimal = selectOptimal(notPlugged)
    cleanupPlates(notPlugged, optimal)

    return optimal


def selectOptimal(plates, **kwargs):
    """Returns the optimal plate to observe."""

    # First tries to select the complete plate with the fewer number
    # of exposures
    completedPlates = [plate for plate in plates if plate.isComplete]
    if len(completedPlates) > 0:
        nExposures = [plate.getTotoroExposures() / plate.priority
                      for plate in completedPlates]
        return completedPlates[np.argmin(nExposures)]

    # If no completed plates, calculates the plate completion using only
    # complete sets
    plateCompletion = np.array(
        [plate.getPlateCompletion() for plate in plates])
    maxPlateCompletion = plateCompletion.max()

    # Selects the plates with maximum completion
    platesMaxCompletion = plates[np.where(
        plateCompletion == maxPlateCompletion)]

    if len(platesMaxCompletion) == 1:
        # If only one plate with maximum completion
        return platesMaxCompletion[0]
    else:
        # If several plates with maximum completion, returns that with the
        # maximum overall completion (including incomplete sets).
        plateCompletionIncomplete = [
            plate.getPlateCompletion(includeIncompleteSets=True)
            for plate in platesMaxCompletion]
        return platesMaxCompletion[plateCompletionIncomplete.argmax()]


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
                if plate.isComplete:
                    break

                lst = site.localSiderealTime(jd)
                if (not utils.isPointInInterval(lst, plateLST, wrapAt=24) or
                        plate.getAltitude(lst) > maxAlt):
                    pass
                else:
                    result = plate.addMockExposure(
                        set=None, startTime=jd, expTime=expTime, silent=True)
                    if result is not False:
                        observedFlag = True
                        result._tmp = True
                    else:
                        observedFlag = False

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


# def getOptimalPlate(plates, JD0, JD1, mode='planner',
#                     prioritiseStarted=True, **kwargs):
#     """Gets the optimal plate to observe between two JDs."""

#     LST0 = site.localSiderialTime(JD0)
#     LST1 = site.localSiderialTime(JD1)

#     visiblePlates = getVisiblePlates(plates, LST0, LST1)

#     if len(visiblePlates) == 0:
#         return None

#     incompletePlates = getIncompletePlates(visiblePlates)

#     if prioritiseStarted:

#         startedPlates = getStartedPlates(incompletePlates)

#         optimumPlateFromStarted = getOptimalPlate(
#             startedPlates, JD0, JD1, mode=mode, prioritiseStarted=False,
#             **kwargs)

#         if optimumPlateFromStarted is not None:
#             return optimumPlateFromStarted

#     simulatedPlates = []

#     for plate in incompletePlates:
#         simulatedPlates.append(simulateExposures(plate, JD0, JD1, **kwargs))

#     platePriority = [plate.getPriority() for plate in simulatedPlates]
#     if 10 in platePriority:
#         optimumPlate = simulatedPlates[platePriority.index(10)]
#         cleanupPlates(simulatedPlates, optimumPlate)
#         return optimumPlate

#     plateCompletion = np.array(
#         [plate.getPlateCompletion() for plate in simulatedPlates])

#     completePlateIdx = np.where(plateCompletion > 1)[0]
#     if len(completePlateIdx) > 0:
#         completedPlates = [simulatedPlates[ii] for ii in completePlateIdx]
#         JD1LastExposure = [plate.getLastExposure().getJD()[1]
#                            for plate in completedPlates]
#         optimumPlate = completedPlates[np.argmin(JD1LastExposure)]
#         cleanupPlates(simulatedPlates, optimumPlate)
#         return optimumPlate

#     platesWithNewExposures = [plate for plate in simulatedPlates
#                               if hasNewExposures(plate)]

#     if len(platesWithNewExposures) == 0:
#         return None

#     plateCompletion = np.array(
#         [plate.getPlateCompletion(includeIncompleteSets=True)
#          for plate in platesWithNewExposures])
#     optimumPlate = platesWithNewExposures[np.argmax(plateCompletion)]
#     cleanupPlates(simulatedPlates, optimumPlate)
#     return optimumPlate

#     # plateNewExposures = np.array(zip(*plateCompletion)[1])
#     # completionTable = table.Table(rows=zip(*plateCompletion)[0],
#     #                               names=['plate', 'completion', 'nSets',
#     #                                      'blueSN2', 'redSN2', 'nNewExp'])

#     # completionTable = completionTable[completionTable['nNewExp'] > 0]

#     # if mode == 'planner':
#     #     priorities = [5] * len(completionTable)

#     # elif mode == 'plugger':
#     #     priorities = []
#     #     for plate in completionTable['plate']:
#     #         platePriority = plate.getPriority()
#     #         priorities.append(platePriority)

#     # completionTable.add_column(
#     #     table.Column(priorities, 'priority'))

#     # if len(completionTable) == 0:
#     #     return (None, [])

#     # completionTable = completionTable[completionTable['priority'] > 1]

#     # if len(completionTable[completionTable['priority'] == 10]):
#     #     completionTable = completionTable[completionTable['priority'] == 10]

#     # if len(completionTable[completionTable['completion'] >= 1.]) > 0:

#     #     idx = np.where(completionTable['completion'] >= 1.)
#     #     selectedStatus = completionTable[idx]
#     #     selectedStatus.sort('nSets')

#     #     if len(selectedStatus[selectedStatus['nSets'] ==
#     #            selectedStatus['nSets'][0]]) > 0:

#     #         idx = np.where(selectedStatus[selectedStatus['nSets'] ==
#     #                        selectedStatus['nSets'][0]])

#     #         selectedStatus = selectedStatus[idx]

#     #         selectedStatus.sort('completion')
#     #         selectedStatus.reverse()

#     #         plate = selectedStatus['plate'][0]
#     #         newExposures = plateNewExposures[
#     #             plateNewExposures[:, 0] == plate][0][1]
#     #         return (plate, newExposures)

#     #     else:

#     #         plate = selectedStatus['plate'][0]
#     #         newExposures = plateNewExposures[
#     #             plateNewExposures[:, 0] == plate][0][1]
#     #         return (plate, newExposures)

#     # else:

#     #     completionTable.sort('completion')

#     #     plate = completionTable['plate'][-1]
#     #     newExposures = plateNewExposures[
#     #         plateNewExposures[:, 0] == plate][0][1]
#     #     return (plate, newExposures)


# def hasNewExposures(plate):
#     """Returns True if a plate has new exposures."""

#     for set in plate.sets:
#         for exp in set.totoroExposures:
#             if exp._tmp:
#                 return True
#     return False


# def getStartedPlates(plates):
#     """Gets plates that are not complete but with completion > 0."""

#     from sdss.internal.manga.Totoro import dbclasses

#     return dbclasses.Plates(
#         [plate for plate in plates
#          if not plate.isComplete and
#             plate.getPlateCompletion(includeIncompleteSets=True) > 0.])


# def getIncompletePlates(plates):
#     """Gets incomplete plates."""

#     from sdss.internal.manga.Totoro import dbclasses

#     incompletePlates = dbclasses.plate.Plates(
#         [plate for plate in plates if not plate.isComplete])

#     return incompletePlates


# def getVisiblePlates(plates, LST0, LST1):
#     """Gets plates that are visible in a range of LST."""

#     from sdss.internal.manga.Totoro import dbclasses

#     return dbclasses.plate.Plates(
#         [plate for plate in plates
#          if utils.getIntervalIntersectionLength(
#              np.array([LST0, LST1]) % 24.,
#              plate.getLSTRange(), wrapAt=24.) > 0])





# def simulateExposures(plate, JD0, JD1, expTime=None, **kwargs):

#     LST0 = site.localSiderialTime(JD0)

#     plateLSTRange = plate.getLSTRange()

#     if not utils.isPointInInterval(LST0, plateLSTRange):
#         return plate

#     currentJD = JD0

#     if expTime is None:
#         expTime = config['exposure']['exposureTime']

#     while True:

#         if currentJD > JD1:
#             return plate

#         currentLST = site.localSiderealTime(currentJD)

#         endJD = currentJD + expTime / config['planner']['efficiency'] / 86400.
#         endLST = site.localSiderealTime(endJD)

#         if not utils.isPointInInterval(
#                 endLST, plate.getLSTRange(), wrapAt=24):
#             return plate

#         maxAlt = config['exposure']['maxAltitude']
#         if plate.getAltitude(currentLST) > maxAlt or \
#                 plate.getAltitude(endLST) > maxAlt:
#             return plate

#         if plate.isComplete:
#             return plate

#         newExposure = plate.addMockExposure(set=None, startTime=currentJD,
#                                             expTime=expTime, silent=True)
#         if newExposure is False:
#             raise TotoroError('something went wrong while creating mock a '
#                               'exposure')

#         newExposure._tmp = True

#         currentJD = endJD
