#!/usr/bin/env python
# encoding: utf-8
"""
scheduler.py

Created by José Sánchez-Gallego on 17 Feb 2014.
Licensed under a 3-clause BSD license.

Revision history:
    17 Feb 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from ..exceptions import TotoroError
from .observingPlan import ObservingPlan
from astropy import time
from .fields import Fields
from astropy.coordinates.angles import Longitude
from astropy.units import degree, hour
from ..core.defaults import *
from .. import log


LONG = Longitude(LONGITUDE(), unit=degree)
ONE_EXPOSURE = EXPTIME() / EFFICIENCY()


class BaseScheduler(object):
    """The base class for the autoscheduler.

    This class provides the common base to plan observations
    and is inherited by Scheduler, Plugger and Planner.

    Parameters
    ----------
    startTime : `astropy.time.Time` or None
        The time at which the scheduling will start. If None,
        the current time is assumed.
    endTime : `astropy.time.Time` or None
        The ending time for the scheduling. If None, the value depends
        on `mode`. If mode='observer', the end of the night will be used;
        if mode='plugger', the end of the next nigh; and, if mode='planner',
        the final date for the survey will be used.
    kwargs
        Any parameter that can be passed to `ObservingPlan`

    """

    def __init__(self, startTime=None, endTime=None,
                 useObservingPlan=True, **kwargs):

        self.longitude = LONG
        self.kwargs = kwargs

        log.debug('Loading observing plan')
        self.observingPlan = ObservingPlan(**self.kwargs)
        self.setStartEndTime(startTime, endTime, **kwargs)

    def setStartEndTime(self, startTime, endTime,
                        scope='survey', **kwargs):
        """Sets the start and end time if they haven't been defined."""

        if startTime is None:
            # Uses the current LST
            startTime = time.Time.now()

        elif isinstance(startTime, time.Time):
            pass

        else:
            TotoroError('startTime must be an astropy.time.Time object.')

        if endTime is None:
            if scope == 'survey':
                endTime = self.observingPlan.getSurveyEnd()
            else:
                raise TotoroError('not yet implemented.')

        elif isinstance(endTime, time.Time):
            pass

        else:
            raise TotoroError('endTime must be an astropy.time.Time object.')

        startTime.lon = self.longitude
        endTime.lon = self.longitude

        # If scope is survey, we only schedule complete days
        if scope == 'survey':

            startTime = self.observingPlan.getClosest(startTime)[0]
            endTime = self.observingPlan.getClosest(endTime)[1]

            self.startTime = time.Time(startTime, format='jd',
                                       scale='utc',
                                       lon=self.longitude)
            self.endTime = time.Time(endTime, format='jd',
                                     scale='utc',
                                     lon=self.longitude)

    def getObservingBlocks(self):
        """Return a table with the observing blocks."""

        return self.observingPlan.getObservingBlocks(
            self.startTime, self.endTime)


class Planner(BaseScheduler):

    def __init__(self, **kwargs):

        self._kwargs = kwargs
        self.reset()

    def getFields(self, **kwargs):
        """Gets a table with the fields that can be scheduled."""
        self.fields = Fields(**kwargs)

        log.debug('Removing done fields')
        self.fields.removeDone()

    def simulate(self):
        """Applies the scheduling logic and returns the planner schedule."""

        for row in self.observingPlan[0:2]:
            startJD = time.Time(row[SURVEY() + '_0'], scale='utc', format='jd')
            endJD = time.Time(row[SURVEY() + '_1'], scale='utc', format='jd')
            self._observeNight(startJD, endJD)

    def _observeNight(self, startJD, endJD):

        def remainingTime(tt):
            return (endJD - tt).sec

        self.fields.plugFields(startJD, endJD)
        currentTime = startJD

        while remainingTime(currentTime) > ONE_EXPOSURE:
            print(currentTime)
            fieldToObserve = self.fields.getOptimumField(currentTime)
            currentTime = fieldToObserve.observe(currentTime, until=endJD)
            print(currentTime)
            if currentTime is None:
                print('Help!!!')

    def reset(self):
        """Reset the class to its original state."""

        super(Planner, self).__init__(**self._kwargs)
        self.getFields(**self._kwargs)
