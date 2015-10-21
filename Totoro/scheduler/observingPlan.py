#!/usr/bin/env python
# encoding: utf-8
"""
observingPlan.py

Created by José Sánchez-Gallego on 11 Dec 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

from astropy.table import Table
from ..core import ConfigObject
from ..exceptions import TotoroError
import os
from astropy import time
from astropy.io import ascii
import warnings


DEFAULT_PLAN_FILE = ConfigObject('defaultPlanFile',
                                 os.path.join(
                                     os.path.dirname(__file__),
                                     'nightly.D.txt'),
                                 'The file with the observing plan')

DEFAULT_PLAN_FORMAT = ConfigObject('defaultPlanFormat',
                                   os.path.join(
                                       os.path.dirname(__file__),
                                       'autoscheduler'),
                                   'The format of the default plan')

OPTIMISED_PLAN_PATTERN = ConfigObject('optimisedPlanPattern',
                                      '{0}.jd',
                                      'The format for the optimised plan')

MJD_COLNAMES = ['MJD', 'APOGEE_0', 'APOGEE_1', 'MaNGA_0', 'MaNGA_1',
                'eBOSS_0', 'eBOSS_1']
AUTOSCHEDULER_COLNAMES = ['col{0}'.format(ii) for ii in range(1, 14)]
AUTOSCHEDULER_VALID_COLNAMES = ['col1', 'col4', 'col5', 'col8',
                                'col9', 'col12', 'col13']


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
                 survey=None, **kwargs):

        self.useOptimisedPlan = useOptimisedPlan
        self.survey = survey
        self.file = None

        if isinstance(plan, Table):
            self._initFromTable(plan)

        if plan is None:
            self.file = DEFAULT_PLAN_FILE()
            self.format = 'autoscheduler'
            self._initFromFile()

        elif isinstance(plan, basestring):
            self.file = plan
            self.format = format.lower()
            self._initFromFile()

        if saveOptimisedPlan and self.file is not None:
            optimisedPlanPath = self.getOptimisedPlanPath(self.file)
            if optimisedPlanPath is not None:
                self.save(optimisedPlanPath)

        if self.survey is not None:
            self.setSurvey(self.survey)

    def _initFromTable(self, plan):
        """Initialises an observing plan from an ObservingPlan instance."""

        self.format = 'jd'
        self.table = plan
        del plan

    def _initFromFile(self):
        """Initialises an observing plan from a file."""

        if self.file is None:
            raise TotoroError('no file to load.')
        elif not os.path.exists(self.file):
            raise TotoroError('file {0} cannot be found.'.format(self.file))

        if self.useOptimisedPlan:
            optimisedPlanPath = self.getOptimisedPlanPath(self.file)
            if os.path.exists(optimisedPlanPath):
                self.file = optimisedPlanPath
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

        self.table = tab
        del tab

    def setSurvey(self, survey):
        """Creates a new instance of the object with only the selected
        survey columns."""

        self.table.keep_columns(['MJD', survey + '_0', survey + '_1'])

    @staticmethod
    def getOptimisedPlanPath(file):
        if not isinstance(file, basestring):
            return None
        configPath = OPTIMISED_PLAN_PATTERN.configPath
        return os.path.join(
            configPath, OPTIMISED_PLAN_PATTERN().format(
                os.path.basename(file)))

    @staticmethod
    def convertSiderealToJD(tab):
        """Converts all the times to JD.

        This function is not very efficient and takes a long
        time to convert all the times for a long observing plan.
        It is recommended to use an input file with the observing
        times already in UTC or JD format.

        """

        from ..utils import lmst2utc

        for row in tab:
            mjd = row[0]
            for nn in range(1, len(row)):
                sidTime = row[nn]
                if sidTime < 0.:
                    continue
                utc = lmst2utc(sidTime, mjd)
                jd = utc.jd
                row[nn] = jd

        return tab

    @staticmethod
    def convertUTCToJD(tab):
        raise TotoroError('convertUTCToJD not yet implemented.')

    def save(self, file):
        """Saves the observing plan to a file."""

        ascii.write(self.table, file, format='fixed_width',
                    delimiter=' ', delimiter_pad=' ')

    def getSurveyEnd(self, survey):
        """Gets the end of survey date."""

        validDates = self[self[survey + '_1'] > 0.0]

        endTime = validDates[-1][survey + '_1']

        tt = time.Time(endTime, format='jd', scale='utc')

        return tt

    def getSurvey(self, survey):
        if survey is None:
            if self.survey is None:
                raise TotoroError('missing survey definition.')
            else:
                return self.survey
        else:
            return survey

    def getStartEndForMJD(self, mjd, survey=None):

        survey = self.getSurvey(survey)
        data = self[self['MJD'] == int(mjd)]

        if len(data) == 0:
            return None

        return data[survey + '_0', survey + '_1'][0]

    def getObservingBlocks(self, startTime, endTime, survey=None):
        """Returns an astropy table with the observation times
        for each night between startTime and endTime."""

        survey = self.getSurvey(survey)

        validDates = self[(self['MJD'] >= int(startTime.mjd)) &
                          (self['MJD'] <= int(endTime.mjd)) &
                          (self[survey + '_0'] > 0.0) &
                          (self[survey + '_1'] > 0.0)]

        if validDates[survey + '_1'][0] < startTime.jd:
            validDates = validDates[1:]
        else:
            validDates[survey + '_0'][0] = startTime.jd

        if validDates[survey + '_0'][-1] > endTime.jd:
            validDates = validDates[0:-1]
        else:
            validDates[survey + '_1'][-1] = endTime.jd

        return validDates[survey + '_0', survey + '_1']

    def keep_columns(self, columns):
        return self.table.keep_columns(columns)

    def __getitem__(self, slice):
        return self.table.__getitem__(slice)

    def __repr__(self):
        return self.table.__repr__()

    def __str__(self):
        return self.table.__str__()
