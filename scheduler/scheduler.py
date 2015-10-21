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
from sdss.internal.manga.Totoro import log, site
from observingPlan import ObservingPlan
from plugger import PluggerScheduler
from planner import PlannerScheduler
from astropy import time


class BaseScheduler(object):
    """The base class for the autoscheduler.

    This class provides the common base to plan observations
    and is inherited by Scheduler, Plugger and Planner.

    """

    def __init__(self, startDate=None, endDate=None, **kwargs):

        self._observingPlan = ObservingPlan()

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

        log.info('Start date: {0}'.format(self.startDate))
        log.info('End date: {0}'.format(self.endDate))


class Planner(BaseScheduler):

    def __init__(self, startDate=None, endDate=None, **kwargs):

        log.info('entering PLANNER mode.')

        super(Planner, self).__init__(startDate=startDate,
                                      endDate=endDate, scope='planner',
                                      **kwargs)

        self._plannerScheduler = PlannerScheduler(self.observingBlocks)
        self._plannerScheduler.schedule(**kwargs)

    #     self.fields = self.getFields(**kwargs)

    # def getFields(self, rejectDrilled=True, **kwargs):
    #     """Gets a table with the fields that can be scheduled."""

    #     from sdss.internal.manga.Totoro import dbclasses

    #     log.info('finding fields with rejectDrilled={0}'
    #              .format(rejectDrilled))
    #     fields = dbclasses.Fields(rejectDrilled=rejectDrilled, **kwargs)

    #     return fields

    # def getExposures(self):
    #     """Returns a table with the simulated exposures."""

    #     if not hasattr(self, 'timelines'):
    #         raise TotoroError('scheduleTimelines must be run before '
    #                           'getExposures')

    #     tableExposures = table.Table(None,
    #                                  names=['manga_tileid', 'RA',
    #                                         'Dec', 'JD0', 'JD1'],
    #                                  dtype=[int, float, float, float, float])

    #     for field in self.fields:
    #         for exposure in field.getValidExposures():
    #             if exposure._tmp is True:
    #                 continue
    #             JD0, JD1 = exposure.getJD()
    #             tableExposures.add_row((field.manga_tileid, field.ra,
    #                                     field.dec, JD0, JD1))

    #     tableExposures.sort('JD0')

    #     return tableExposures

    # def scheduleTimelines(self):
    #     """Creates simulated timelines for the observing blocks."""

    #     self.timelines = Timelines(
    #         self.observingBlocks, plates=self.fields, mode='planner')
    #     self.timelines.schedule()


class Plugger(object):

    def __init__(self, startDate, endDate, **kwargs):

        # Not necessary for now.
        # super(Plugger, self).__init__(startDate=startDate, endDate=endDate,
        #                               scope='plugger', **kwargs)
        # if len(self.observingBlocks) == 0:
        #     raise TotoroError('no observing blocks found.')
        # elif len(self.observingBlocks) > 1:
        #     raise TotoroError('Plugger must be run only for one night. '
        #                       'Multiple observing blocks found.')

        assert startDate < endDate

        self.startDate = startDate
        self.endDate = endDate
        log.info('Start date: {0}'.format(self.startDate))
        log.info('End date: {0}'.format(self.endDate))
        log.info('Scheduling {0:.2f} hours'.format(
                 (self.endDate - self.startDate)*24.))

        self.plates = self.getPlatesAtAPO()

    def getPlatesAtAPO(self, rejectComplete=False):

        from sdss.internal.manga.Totoro import dbclasses

        log.info('getting plates at APO with rejectComplete={0}'.format(
                 rejectComplete))

        plates = dbclasses.getAtAPO(onlyIncomplete=rejectComplete)

        log.info('plates found: {0}'.format(len(plates)))

        return plates

    def getOutput(self, **kwargs):

        pluggerSchedule = PluggerScheduler(
            self.plates, self.startDate, self.endDate, **kwargs)

        return pluggerSchedule.carts


class Nightly(object):

    def __init__(self, startDate=None, endDate=None, plates=None, **kwargs):

        log.info('entering NIGHTLY mode.')

        # Not necessary, as we don't really use the schedule for now.
        # super(Nightly, self).__init__(startDate=startDate,
        #                               endDate=endDate, scope='nightly',
        #                               **kwargs)

        self.startDate = startDate
        self.endDate = endDate
        self.currentDate = time.Time.now().jd

        if plates is None:
            self.plates = self.getPlates(**kwargs)
        else:
            log.info('using programmatic input for plates.')
            log.info('{0} input plates'.format(len(plates)))
            self.plates = plates

    def getPlates(self, **kwargs):
        """Gets the plugged plates."""

        from sdss.internal.manga.Totoro import dbclasses

        plates = dbclasses.getPlugged(**kwargs)

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
