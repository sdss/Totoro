#!/usr/bin/env python
# encoding: utf-8
"""
observingPlan.py

Created by José Sánchez-Gallego on 11 Dec 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

from __future__ import division
from __future__ import print_function
from astropy.table import Table, Column
from ..exceptions import TotoroError
from astropy import time
from astropy.io import ascii
from ..core.defaults import SURVEY, DEFAULT_PLAN_FILE, \
    OPTIMISED_PLAN_PATTERN, MJD_COLNAMES, AUTOSCHEDULER_VALID_COLNAMES
import numpy as np
import os


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

    def __init__(self, plan=None, format='jd',
                 useOptimisedPlan=True, saveOptimisedPlan=True,
                 **kwargs):

        self._useOptimisedPlan = useOptimisedPlan
        self.file = None
        self.isFileOptimised = False
        self.data = None

        if 'survey' in kwargs:
            self.survey = kwargs['survey']
        else:
            self.survey = SURVEY()

        if isinstance(plan, Table):
            self._initFromTable(plan)

        elif plan is None:
            self.file = DEFAULT_PLAN_FILE()
            self.format = 'autoscheduler'
            self._initFromFile()

        elif isinstance(plan, basestring):
            self.file = plan
            self.format = format.lower()
            self._initFromFile()

        if saveOptimisedPlan and self.file is not None:
            optimisedPlanPath = self.getOptimisedPlanPath()
            if optimisedPlanPath is not None:
                self.save(optimisedPlanPath)

        self.data = self.data['MJD', self.survey + '_0', self.survey + '_1']
        self.addRunDayCol()
        self.data = self.data[self.data[self.survey + '_0'] != -1]

    def _initFromTable(self, plan):
        """Initialises an observing plan from an ObservingPlan instance."""

        self.format = 'jd'
        self.data = plan

    def _initFromFile(self):
        """Initialises an observing plan from a file."""

        if self.file is None:
            raise TotoroError('no file to load.')
        elif not os.path.exists(self.file):
            raise TotoroError('file {0} cannot be found.'.format(self.file))

        if self._useOptimisedPlan:
            optimisedPlanPath = self.getOptimisedPlanPath()
            if os.path.exists(optimisedPlanPath):
                self.file = optimisedPlanPath
                self.isFileOptimised = True
                self.format = 'jd'

        if 'autoscheduler' in self.format:
            tab = ascii.read(self.file, format='no_header', data_start=1)
            tab = tab[AUTOSCHEDULER_VALID_COLNAMES]
            for nn, col in enumerate(tab.colnames):
                tab.rename_column(col, MJD_COLNAMES[nn])
        else:
            tab = ascii.read(self.file, format='no_header', data_start=1,
                             names=MJD_COLNAMES)

        if self.format in ['autoscheduler', 'sidereal']:
            tab = self.convertSiderealToJD(tab)
        elif self.format in ['autoscheduler.utc', 'utc']:
            tab = self.convertUTCToJD(tab)
        elif self.format == 'autoscheduler.jd':
            pass
        else:
            TotoroError('format not understood.')

        self.data = tab

    def getOptimisedPlanPath(self):
        if self.isFileOptimised:
            return self.file
        elif not isinstance(self.file, basestring):
            return None
        else:
            configPath = OPTIMISED_PLAN_PATTERN.configPath
            return os.path.join(
                configPath, OPTIMISED_PLAN_PATTERN().format(
                    os.path.basename(self.file)))

    @staticmethod
    def convertSiderealToJD(tab):
        """Converts all the times to JD.

        This function is not very efficient and takes a long
        time to convert all the times for a long observing plan.
        It is recommended to use an input file with the observing
        times already in UTC or JD format.

        """

        from ..utils import lmst2time

        for row in tab:
            mjd = row[0]
            for nn in range(1, len(row)):
                sidTime = row[nn]
                if sidTime < 0.:
                    continue
                utc = lmst2time(sidTime, mjd)
                jd = utc.jd
                row[nn] = jd

        return tab

    @staticmethod
    def convertUTCToJD(tab):
        raise TotoroError('convertUTCToJD not yet implemented.')

    def save(self, file):
        """Saves the observing plan to a file."""

        ascii.write(self.data, file, format='fixed_width',
                    delimiter=' ', delimiter_pad=' ')

    def addRunDayCol(self):
        """Adds a column with the night within the run."""

        nDay = 1
        ll = []

        for row in self.data:
            if row[self.survey + '_0'] != -1:
                ll.append(nDay)
                nDay += 1
            else:
                ll.append(-1)
                nDay = 1

        self.data.add_column(Column(data=ll, name='RUN_DAY', dtype=int))

    def getSurveyEnd(self):
        """Gets the end of survey date."""

        survey = self.survey
        validDates = self.data[self.data[survey + '_1'] > 0.0]
        endTime = validDates[-1][survey + '_1']
        tt = time.Time(endTime, format='jd', scale='utc')
        return tt

    def getClosest(self, dd, survey=SURVEY()):

        tt = self.data[self.data[self.survey + '_1'] >= dd.jd]

        idxStart = (np.abs(tt[self.survey + '_0'] - dd.jd)).argmin()
        idxEnd = (np.abs(tt[self.survey + '_1'] - dd.jd)).argmin()

        return (tt[idxStart][self.survey + '_0'],
                tt[idxEnd][self.survey + '_1'])

    # def getObservingBlocks(self, startTime, endTime, survey=None):
    #     """Returns an astropy table with the observation times
    #     for each night between startTime and endTime."""

    #     survey = self.getSurvey(survey)

    #     validDates = self[(self['MJD'] >= int(startTime.mjd)) &
    #                       (self['MJD'] <= int(endTime.mjd)) &
    #                       (self[survey + '_0'] > 0.0) &
    #                       (self[survey + '_1'] > 0.0)]

    #     if validDates[survey + '_1'][0] < startTime.jd:
    #         validDates = validDates[1:]
    #     else:
    #         validDates[survey + '_0'][0] = startTime.jd

    #     if validDates[survey + '_0'][-1] > endTime.jd:
    #         validDates = validDates[0:-1]
    #     else:
    #         validDates[survey + '_1'][-1] = endTime.jd

    #     return validDates[survey + '_0', survey + '_1']

    # def __getitem__(self, slice):
    #     return Table.__getitem__(self, slice)

    def __repr__(self):
        return self.data.__repr__()

    def __str__(self):
        return self.data.__str__()

    def __getitem__(self, slice):
        return self.data[slice]
