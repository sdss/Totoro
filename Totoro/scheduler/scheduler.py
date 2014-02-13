#!/usr/bin/env python
# encoding: utf-8
"""
Scheduler.py

Created by José Sánchez-Gallego on 11 Dec 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

from ..exceptions import TotoroError
from observingPlan import ObservingPlan
from astropy import time
from ..core import ConfigObject
from astropy.coordinates.angles import Longitude
from astropy.units import degree

APO_LONGITUDE = 254.179722
LONGITUDE = ConfigObject('longitude', APO_LONGITUDE,
                         'The longitude of the observatory')

SURVEY = ConfigObject('survey', 'MaNGA',
                      'The survey to schedule')


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

        self.longitude = Longitude(LONGITUDE(), unit=degree)
        self.kwargs = kwargs

        if 'survey' in self.kwargs:
            self.survey = self.kwargs['survey']
        else:
            self.survey = SURVEY()
            self.kwargs['survey'] = SURVEY()

        self.startTime = startTime
        self.endTime = endTime

        self.observingPlan = ObservingPlan(**kwargs)

    def setStartEndTime(self, scope='night'):
        """Sets the start and end time if they haven't been defined."""

        if self.startTime is None:
            # Uses the current LST
            self.startTime = time.Time.now()

        elif isinstance(self.startTime, time.Time):
            pass

        else:
            TotoroError('startTime must be an astropy.time.Time object.')

        self.startTime.lon = self.longitude

        if self.endTime is None:
            if scope == 'survey':
                self.endTime = self.observingPlan.getSurveyEnd(self.survey)
            else:
                raise TotoroError('not yet implemented.')

        elif isinstance(self.endTime, time.Time):
            pass

        else:
            raise TotoroError('endTime must be an astropy.time.Time object.')

        self.endTime.lon = self.longitude

        # If scope is survey, we only schedule complete days
        if scope == 'survey':
            startTimeJD = self.observingPlan.getStartEndForMJD(
                self.startTime.mjd, self.survey)[0]

            endTimeJD = self.observingPlan.getStartEndForMJD(
                self.endTime.mjd, self.survey)[1]

            self.startTime = time.Time(startTimeJD, format='jd',
                                       scale='utc',
                                       lon=self.longitude)
            self.endTime = time.Time(endTimeJD, format='jd',
                                     scale='utc',
                                     lon=self.longitude)

    def getObservingBlocks(self):
        """Return a table with the observing blocks."""

        return self.observingPlan.getObservingBlocks(
            self.startTime, self.endTime, self.survey)


class Planner(BaseScheduler):

    def __init__(self, **kwargs):

        super(Planner, self).__init__(**kwargs)
        self.setStartEndTime(scope='survey')
        self.observingBlocks = self.getObservingBlocks()

    def getTiles(self):
        """Gets a table with the tiles/plates that can be scheduled."""




