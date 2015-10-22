#!/usr/bin/env python
# encoding: utf-8
"""
intervals.py

Created by José Sánchez-Gallego on 4 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    4 Aug 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import numpy as np
from Totoro.exceptions import TotoroError
import itertools


def getIntervalIntersectionLength(aa, bb, wrapAt=360):
    """Returns the length of the instersection between two intervals aa and bb.
    """

    intersection = getIntervalIntersection(aa, bb, wrapAt=wrapAt)

    if intersection is False:
        return 0.0
    else:
        if wrapAt is None:
            return (intersection[1] - intersection[0])
        else:
            return (intersection[1] - intersection[0]) % wrapAt


def getIntervalIntersection(aa, bb, wrapAt=360):
    """Returns the intersection between two intervals."""

    if wrapAt is None:
        if bb[1] - bb[0] > aa[1] - aa[0]:
            aa, bb = bb, aa
    else:
        if (bb[1] - bb[0]) % wrapAt > (aa[1] - aa[0]) % wrapAt:
            aa, bb = bb, aa

    if (isPointInInterval(bb[0], aa, wrapAt=wrapAt) and
            isPointInInterval(bb[1], aa, wrapAt=wrapAt)):
        return np.array([bb[0], bb[1]])

    if (not isPointInInterval(bb[0], aa, wrapAt=wrapAt) and
            not isPointInInterval(bb[1], aa, wrapAt=wrapAt)):
        return False

    if isPointInInterval(bb[0], aa, wrapAt=wrapAt):
        return np.array([bb[0], aa[1]])

    if isPointInInterval(bb[1], aa, wrapAt=wrapAt):
        return np.array([aa[0], bb[1]])


def isPointInInterval(point, ival, wrapAt=360):
    """Returns True if point in interval."""

    if wrapAt is None:
        return point >= ival[0] and point <= ival[1]
    else:
        return (point - ival[0]) % wrapAt <= (ival[1] - ival[0]) % wrapAt


def isIntervalInsideOther(aa, bb, wrapAt=360, onlyOne=False):
    """Checks if the interval aa (a numpy.ndarray of length 2) is inside bb."""

    if wrapAt is None:
        p1 = (aa[0] >= bb[0]) and (aa[0] <= bb[1])
        p2 = (aa[1] >= bb[0]) and (aa[1] <= bb[1])
    else:
        p1 = ((aa[0] - bb[0]) % wrapAt < (bb[1]-bb[0]) % wrapAt)
        p2 = ((aa[1] - bb[0]) % wrapAt < (bb[1]-bb[0]) % wrapAt)

    if p1 and p2:
        return True
    elif onlyOne and (p1 or p2):
        return True

    return False


def intervalLength(aa, wrapAt=360.):
    if wrapAt is None:
        return (aa[1] - aa[0])
    else:
        return (aa[1] - aa[0]) % wrapAt


def getMinMaxIntervalSequence(intervals, wrapAt=360):

    if intervals.shape[0] < 1:
        raise TotoroError('input needs to be an Nx2 array with N>=1.')
    elif intervals.shape[0] == 1:
        return np.array(intervals[0])

    firstTwo = intervals[0:2, :]
    rest = intervals[2:, :]

    interval1 = firstTwo[0, :]
    length1 = intervalLength(interval1, wrapAt=wrapAt)

    interval2 = firstTwo[1, :]
    length2 = intervalLength(interval2, wrapAt=wrapAt)

    if (isIntervalInsideOther(interval1, interval2) or
            isIntervalInsideOther(interval2, interval1)):
        newInterval = interval1 if length1 > length2 else interval2
    else:
        tmpInt1 = [interval1[0], interval2[1]]
        tmpLength1 = intervalLength(tmpInt1, wrapAt=wrapAt)
        tmpInt2 = [interval2[0], interval1[1]]
        tmpLength2 = intervalLength(tmpInt2, wrapAt=wrapAt)

        if (isPointInInterval(interval1[0], interval2) or
                isPointInInterval(interval1[1], interval2)):
            newInterval = tmpInt1 if tmpLength1 > tmpLength2 else tmpInt2
        else:
            newInterval = tmpInt1 if tmpLength1 < tmpLength2 else tmpInt2

    if rest.shape[0] == 0:
        return newInterval
    else:
        return getMinMaxIntervalSequence(
            np.append([newInterval], rest, axis=0))


def calculateMean(interval, wrapAt=360.):
    if wrapAt is None:
        return (interval[0] + (interval[1] - interval[0]) / 2.)
    else:
        return ((interval[0] + ((interval[1] - interval[0]) % wrapAt) / 2.) %
                wrapAt)


def getIntervalFromPoints(points, wrapAt=360.):

    if wrapAt is None:
        points = np.array(points)
    else:
        points = np.array(points) % wrapAt

    if len(points) == 1:
        return points

    validExtremes = []
    for permutation in itertools.permutations(points, 2):
        isValid = True
        for point in points:
            if not isPointInInterval(point, permutation):
                isValid = False
        if isValid:
            validExtremes.append(permutation)

    lengths = np.array(
        [intervalLength(interval) for interval in validExtremes])

    return validExtremes[np.argmin(lengths)]


def removeInterval(master, interToRemove, wrapAt=360.):
    """Removes an interval within another interval or set of intervals"""

    master = np.atleast_2d(master)

    if master.size == 0:
        return master

    if master.shape[0] == 1:
        if isIntervalInsideOther(master[0], interToRemove, wrapAt=wrapAt):
            return np.array([])

    newMaster = []
    for interval in master:

        intersection = getIntervalIntersection(interval, interToRemove,
                                               wrapAt=wrapAt)

        if intersection is False:
            newMaster.append(interval)
            continue

        if intersection[0] == interval[0] and intersection[1] == interval[1]:
            pass
        elif intersection[0] == interval[0]:
            newMaster.append([intersection[1], interval[1]])
        elif intersection[1] == interval[1]:
            newMaster.append([interval[0], intersection[0]])
        else:
            newMaster.append([interval[0], intersection[0]])
            newMaster.append([intersection[1], interval[1]])

    if wrapAt is None:
        newMaster = np.array(newMaster)
    else:
        newMaster = np.array(newMaster) % wrapAt

    if newMaster.shape[0] == 1:
        newMaster = newMaster[0]

    return newMaster


def splitInterval(interval, wrapAt=360.):
    """Splits an interval in a list of intervals that do not
    wrap around wrapAt. For instance, the interval [300, 120] would be
    returned as [[300, 360], [0, 120]]."""

    assert len(interval) == 2, 'interval must have length 2'

    interval = np.array(interval) % wrapAt

    if interval[0] < interval[1]:
        return interval

    return np.array([[interval[0], wrapAt], [0, interval[1]]])
