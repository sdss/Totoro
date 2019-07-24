#!/usr/bin/env python
# encoding: utf-8
"""
observingPlan.py

Created by José Sánchez-Gallego on 11 Dec 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

from __future__ import division, print_function

import os

import numpy as np
from astropy import table, time
from six import string_types

from Totoro import config, exceptions, log, readPath


def getScheduleFile():

    try:
        # Gets the master scheduler schedule file
        import autoscheduler
        autoschedulerDir = os.path.dirname(autoscheduler.__file__)
        schedule = os.path.join(autoschedulerDir, '../../schedules', 'Sch_APO_base.dat')
    except Exception:
        # If the autoscheduler is not available, uses the local schedule file
        schedule = readPath(config['observingPlan']['fallBackSchedule'])
        log.warning(
            'The master autoscheduler could not be found. '
            'Using a local copy of the schedule. BE CAREFUL! '
            'THIS FILE MIGHT BE OUTDATED!', exceptions.TotoroUserWarning)

    return schedule


class ObservingPlan(object):
    """The survey-wide observing plan object.

    Parameters
    ----------
    plan : str or `astropy.tableTable`
        A file pth containing the observing plan or an astropy table
        with the data.
    format : str
        The input format of the plan. The options are 'autoscheduler',
        'autoscheduler.utc', 'autoscheduler.jd', 'sidereal', 'utc' or 'jd'.
    useOptimisedPlan : bool
        It True, tries to read the optimised file (in JD format) for the
        input plan. This notably improves performance.
    saveOptimisedPlan : bool
        If True, once the input observing plan has been converted to JD,
        saves it to the configuration directory. This optimised file will
        later be read if useOptimisedPlan is True.
    survey : str
        If None, the observing plan for all plans will be kept, but survey
        will need to be defined for some methods.

    """

    def __init__(self, schedule=config['observingPlan']['schedule'], **kwargs):

        if isinstance(schedule, string_types) and schedule.lower() == 'none':
            schedule = getScheduleFile()
        else:
            schedule = readPath(schedule)

        self.scheduleFile = os.path.realpath(schedule)
        scheduleToPrint = self.scheduleFile[:15] + '...' + \
            self.scheduleFile[-45:] if len(self.scheduleFile) > 70 else \
            self.scheduleFile

        if not os.path.exists(self.scheduleFile):
            raise exceptions.TotoroError('schedule {0} not found'.format(scheduleToPrint))

        plan = table.Table.read(self.scheduleFile, format='ascii.no_header')

        self.plan = plan.copy()
        self.plan.keep_columns(['col1', 'col11', 'col12'])
        self.plan.rename_column('col1', 'JD')
        self.plan.rename_column('col11', 'JD0')
        self.plan.rename_column('col12', 'JD1')

        self.plan.add_column(table.Column(np.zeros(len(self.plan)), name='Position', dtype=int))
        for ii, row in enumerate(plan):
            if row['col11'] == 0:
                continue
            if row['col11'] > row['col5'] and row['col11'] > row['col9']:
                self.plan['Position'][ii] = 2
            else:
                self.plan['Position'][ii] = 1

        log.debug('observing plan {0} loaded.'.format(scheduleToPrint))

        self.addRunDayCol()
        self.plan = self.plan[(self.plan['JD0'] > 0) & (self.plan['JD1'] > 0)]

    def addRunDayCol(self):
        """Adds a column with the night within the run."""

        nDay = 1
        nRun = 1
        ll = []
        run = []

        for row in self.plan:
            if row['JD0'] != 0.:
                ll.append(nDay)
                run.append(nRun)
                nDay += 1
            else:
                ll.append(-1)
                run.append(-1)
                nDay = 1
                nRun += 1

        run = np.array(run)
        for nn, value in enumerate(np.unique(run[run > 0])):
            run[run == value] = nn + 1

        self.plan.add_column(table.Column(data=run, name='RUN', dtype=int))
        self.plan.add_column(table.Column(data=ll, name='RUN_DAY', dtype=int))

    def getSurveyStart(self):
        """Gets the start of survey."""

        validDates = self.plan[self.plan['JD0'] > 0.0]
        startTime = validDates[0]['JD0']
        return startTime

    def getSurveyEnd(self):
        """Gets the end of survey."""

        validDates = self.plan[self.plan['JD1'] > 0.0]
        endTime = validDates[-1]['JD1']
        return endTime

    def getClosest(self, dd):

        tt = self.plan[(self.plan['JD0'] < dd) & (self.plan['JD1'] > dd)]

        if len(tt) == 1:
            return tt

        tt = self.plan[self.plan['JD'] <= dd]

        return tt[-1]

    def getJD(self, jd=None):

        if jd is None:
            jd = time.Time.now().jd

        jd = int(jd)

        if jd not in self.plan['JD']:
            log.warning('JD={0} not found in schedule'.format(jd), exceptions.TotoroUserWarning)
            return (None, None)

        night = self.plan[self.plan['JD'] == jd]

        return (night['JD0'][0], night['JD1'][0])

    def getObservingBlocks(self, startDate, endDate):
        """Returns an astropy table with the observation dates
        for each night between startDate and endDate."""

        validDates = self.plan[(self.plan['JD1'] >= startDate) & (self.plan['JD0'] <= endDate) &
                               (self.plan['JD0'] > 0.0) & (self.plan['JD1'] > 0.0)]

        if (startDate is None and endDate is None) or (len(validDates) == 0):
            log.info('no observing blocks selected.')
            return validDates

        if startDate > validDates['JD0'][0]:
            validDates['JD0'][0] = startDate

        if endDate < validDates['JD1'][-1]:
            validDates['JD1'][-1] = endDate

        totalTime = np.sum([row['JD1'] - row['JD0'] for row in validDates]) * 24

        log.info(('{0} blocks (days) selected, '
                  'making a total of {1:.2f} hours').format(len(validDates), totalTime))

        return validDates

    def _createObservingBlock(self, startDate, endDate):
        """Creates an observing block from start and end dates. Mainly for the
        plugger."""

        tmpPlan = self.plan.copy()
        tmpPlan.add_row((int(startDate), startDate, endDate, -1, -1))

        totalTime = (endDate - startDate) * 24.
        log.info(('1 block (days) selected, ' 'making a total of {0:.2f} hours').format(totalTime))

        return table.Table(tmpPlan[len(tmpPlan) - 1])

    def getRun(self, startDate=None):

        jd = int(startDate if startDate is not None else time.Time.now().jd)

        if jd not in self.plan['JD']:
            raise exceptions.TotoroError('JD={0} not found in schedule'.format(jd))

        firstDay = self.plan[self.plan['JD'] == jd]
        runNumber = firstDay['RUN']
        lastDay = self.plan[self.plan['RUN'] == runNumber][-1]

        return (firstDay['JD0'], lastDay['JD1'])

    def getMJD(self, mjd):

        mjd = int(mjd)

        jd = 2400000. + int(mjd)
        row = self.plan[self.plan['JD'] == jd]

        if len(row) == 1:
            return (row['JD0'][0], row['JD1'][0])
        else:
            return None

    def getPosition(self, jd):
        """Returns the position of the block in which a certain JD falls."""

        try:
            pos = int(self.getClosest(jd)['Position'])
        except Exception:
            log.warning('Cannot found block position for JD={0}'.format(jd),
                        exceptions.TotoroUserWarning)
            pos = None

        return pos

    def __repr__(self):
        return self.plan.__repr__()

    def __str__(self):
        return self.plan.__str__()

    def __getitem__(self, slice):
        return self.plan[slice]
