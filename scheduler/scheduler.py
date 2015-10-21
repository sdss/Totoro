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
from astropysics import obstools
from sdss.internal.manga.Totoro.exceptions import TotoroError
from timeline import Timelines
from sdss.internal.manga.Totoro import log, site
from sdss.internal.manga.Totoro.scheduler import observingPlan
from astropy import time


class BaseScheduler(object):
    """The base class for the autoscheduler.

    This class provides the common base to plan observations
    and is inherited by Scheduler, Plugger and Planner.

    """

    def __init__(self, startDate=None, endDate=None, **kwargs):

        if observingPlan is None:
            raise TotoroError('observing plan not found. Not possible to '
                              'create an instance of BaseScheduler.')

        self._observingPlan = observingPlan
        self.site = site
        self._setStartEndDate(startDate, endDate, **kwargs)

        self.observingBlocks = self._observingPlan.getObservingBlocks(
            self.startDate, self.endDate)

    def _setStartEndDate(self, startDate, endDate, scope='planner', **kwargs):
        """Sets the start and end date if they haven't been defined."""

        if scope == 'planner':
            if startDate is None or endDate is None:
                raise TotoroError('planner requires both start and end dates.')

        elif scope == 'nightly':
            if startDate is None or endDate is None:
                startDate, endDate = self._observingPlan.getJD(startDate)

        elif scope == 'plugger':
            startDate, endDate = self._observingPlan.getJD(
                int(time.Time.now().jd) + 1)

        self.startDate = startDate
        self.endDate = endDate
        self.currentDate = obstools.calendar_to_jd(None)

        log.debug('Start date: {0}'.format(self.startDate))
        log.debug('End date: {0}'.format(self.endDate))


class Planner(BaseScheduler):

    def __init__(self, startDate=None, endDate=None, **kwargs):

        log.info('entering PLANNER mode.')

        super(Planner, self).__init__(startDate=startDate,
                                      endDate=endDate, scope='planner',
                                      **kwargs)

        self.fields = self.getFields(**kwargs)

    def getFields(self, rejectDrilled=True, **kwargs):
        """Gets a table with the fields that can be scheduled."""

        from sdss.internal.manga.Totoro import dbclasses

        log.info('finding fields with rejectDrilled={0}'.format(rejectDrilled))
        fields = dbclasses.Fields(rejectDrilled=rejectDrilled, **kwargs)

        return fields

    def scheduleTimelines(self):
        """Creates simulated timelines for the observing blocks."""

        self.timelines = Timelines(
            self.observingBlocks, plates=self.fields, mode='planner')
        self.timelines.schedule()


class Plugger(BaseScheduler):

    def __init__(self, date=None, **kwargs):

        super(Plugger, self).__init__(startDate=date, endDate=None,
                                      scope='plugger', **kwargs)

        self.plates = self.getPlatesAtAPO()

    def getPlatesAtAPO(self, rejectComplete=True):

        from sdss.internal.manga.Totoro import dbclasses

        log.info('getting plates at APO with rejectComplete={0}'.format(
                 rejectComplete))

        plates = dbclasses.Plates.getAtAPO(onlyIncomplete=True)

        log.info('plates found: {0}'.format(len(plates)))

        return plates

    def getOutput(self):

        self.timelines = Timelines(
            self.observingBlocks, plates=self.plates, mode='plugger')

        self.timelines.schedule()


class Nightly(BaseScheduler):

    def __init__(self, startDate=None, endDate=None, plates=None, **kwargs):

        log.info('entering NIGHTLY mode.')

        super(Nightly, self).__init__(startDate=startDate,
                                      endDate=endDate, scope='nightly',
                                      **kwargs)

        if plates is None:
            self.plates = self.getPlates(**kwargs)
        else:
            log.info('using programmatic input for plates.')
            log.info('{0} input plates'.format(len(plates)))
            self.plates = plates

    def getPlates(self, **kwargs):
        """Gets the plugged plates."""

        from sdss.internal.manga.Totoro import dbclasses

        plates = dbclasses.Plates.getPlugged(**kwargs)

        if len(plates) == 0:
            log.info('no plugged plates found.')
        else:
            log.info('found {0} plugged plates'.format(len(plates)))

        return plates

    def printTabularOutput(self):
        """Prints a series of tables with information about the schedule."""

        from sdss.internal.manga.Totoro import output

        output.printTabularOutput(self.plates)

    def getOutput(self, format='dict'):
        """Returns the nightly output in the selected format."""

        from sdss.internal.manga.Totoro import output

        if format == 'table':
            return output.getTabularOutput(self.plates)
        else:
            return output.getNightlyOutput(self, format=format)
