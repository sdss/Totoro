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
from sdss.internal.manga.Totoro.exceptions import TotoroError
from sdss.internal.manga.Totoro import log, site, config
from sdss.internal.manga.Totoro.scheduler import observingPlan
from sdss.internal.manga.Totoro.utils import intervals
from plugger import PluggerScheduler
from planner import PlannerScheduler
from astropy import time
import numpy as np


__ALL__ = ['BaseScheduler', 'Planner', 'Plugger']


class BaseScheduler(object):
    """The base class for the autoscheduler.

    This class provides the common base to plan observations
    and is inherited by Scheduler, Plugger and Planner.

    """

    def __init__(self, startDate=None, endDate=None, **kwargs):

        self._observingPlan = observingPlan

        if self._observingPlan is None:
            raise TotoroError('observing plan not found. Not possible to '
                              'create an instance of BaseScheduler.')

        self.site = site
        self._setStartEndDate(startDate, endDate, **kwargs)

        scope = kwargs.get('scope', 'planner')

        if scope != 'plugger':
            self.observingBlocks = self._observingPlan.getObservingBlocks(
                self.startDate, self.endDate)
        else:
            self.observingBlocks = self._observingPlan._createObservingBlock(
                self.startDate, self.endDate)

    def _setStartEndDate(self, startDate, endDate, scope='planner', **kwargs):
        """Sets the start and end date if they haven't been defined."""

        if scope == 'planner':
            if startDate is None:
                startDate = time.Time.now().jd
            if startDate is None or endDate is None:
                raise TotoroError('planner end date.')

        elif scope == 'nightly':
            if startDate is None or endDate is None:
                startDate, endDate = self._observingPlan.getJD(startDate)

        elif scope == 'plugger':
            if startDate >= endDate:
                raise TotoroError('startDate must be < endDate')

        self.startDate = startDate
        self.endDate = endDate
        self.currentDate = time.Time.now().jd

        log.info('start date: {0}'.format(self.startDate))
        log.info('end date: {0}'.format(self.endDate))


class Planner(BaseScheduler):

    def __init__(self, startDate=None, endDate=None, **kwargs):

        log.info('entering PLANNER mode.')

        super(Planner, self).__init__(startDate=startDate,
                                      endDate=endDate, scope='planner',
                                      **kwargs)

        self._plannerScheduler = PlannerScheduler(self.observingBlocks,
                                                  **kwargs)
        self._plannerScheduler.schedule(**kwargs)


class Plugger(object):

    def __init__(self, startDate, endDate, **kwargs):

        assert startDate < endDate

        self.startDate = startDate
        self.endDate = endDate
        log.info('Start date: {0}'.format(self.startDate))
        log.info('End date: {0}'.format(self.endDate))
        log.info('Scheduling {0:.2f} hours'.format(
                 (self.endDate - self.startDate)*24.))

        self.plates = self.getPlatesAtAPO(**kwargs)

    def getPlatesAtAPO(self, rejectComplete=False, onlyMarked=False, **kwargs):

        from sdss.internal.manga.Totoro import dbclasses

        onlyVisiblePlates = kwargs.pop('onlyVisiblePlates',
                                       config['plugger']['onlyVisiblePlates'])

        assert isinstance(onlyVisiblePlates, int), \
            'onlyVisiblePlates must be a boolean'

        log.info('getting plates at APO with rejectComplete={0}, '
                 'onlyMarked={1}'.format(rejectComplete, onlyMarked))

        if onlyVisiblePlates:
            lstRange = site.localSiderealTime([self.startDate, self.endDate])
            window = config['plateVisibilityMaxHalfWindowHours']
            raRange = np.array([(lstRange[0] - window) * 15.,
                                (lstRange[1] + window) * 15.])

            log.info('selecting plates with RA in range {0}'
                     .format(str(raRange % 360)))

            raRange = intervals.splitInterval(raRange, 360.)

        else:
            raRange = None

        plates = dbclasses.getAtAPO(onlyIncomplete=rejectComplete,
                                    onlyMarked=onlyMarked,
                                    rejectLowPriority=True,
                                    fullCheck=False, raRange=raRange)

        log.info('plates found: {0}'.format(len(plates)))

        return plates

    def getOutput(self, **kwargs):

        pluggerSchedule = PluggerScheduler(
            self.plates, self.startDate, self.endDate, **kwargs)

        return pluggerSchedule.carts
