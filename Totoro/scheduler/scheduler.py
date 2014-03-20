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
from astropy.units import degree
from ..core.defaults import LONGITUDE, SURVEY, EXPTIME
from collections import OrderedDict
from ..utils.conversion import utc2lmst


LONG = Longitude(LONGITUDE(), unit=degree)


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

        super(Planner, self).__init__(**kwargs)
        self.getFields(**kwargs)
        self.schedule = OrderedDict()

    def getFields(self, **kwargs):
        """Gets a table with the fields that can be scheduled."""
        self.fields = Fields(**kwargs)
        self.fields.removeCompleted()

    def getSchedule(self):
        """Applies the scheduling logic and returns the planner schedule."""

        self.schedule = OrderedDict()
        for row in self.observingPlan:
            mjd = row['MJD']
            startJD = row[SURVEY() + '_0']
            endJD = row[SURVEY() + '_1']
            self.observe(mjd, startJD, endJD)

    def observe(self, mjd, startJD, endJD):

        lstStart = utc2lmst(startJD)
        lstEnd = utc2lmst(endJD)
        currentTime = lstStart

        # Remaining time in minutes
        remainingTime = (lstEnd - lstStart).sec / 60.

    def reset(self):
        """Reset the class to its original state."""
        self.getFields()
