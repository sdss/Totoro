#!/usr/bin/env python
# encoding: utf-8
"""
queries.py

Created by José Sánchez-Gallego on 31 Oct 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

from dataModel import *
from sqlalchemy.orm.attributes import InstrumentedAttribute

__all__ = ['ExposureQuery']


class ExposureQuery(object):
    """Queries the platedb.exposure table.

    Creates a query around the platedb.exposure table, allowing to perform
    some of the most used queries.

    Parameters
    ----------
    query : `query` or None
        If a `query` is specified, all the methods will act on it. Otherwise,
        creates a simple query of the table exposure.

    Examples
    --------

      Examples ::

        qq = ExposureQuery()
        qq.filterScience()

        for row in qq.query[0:10]:
            print row.pk, row.exposureStatus.label, row.survey.label

        2740 Good BOSS
        2741 Good BOSS
        2742 Good BOSS
        2743 Good BOSS
        2744 Good BOSS
        2714 Good BOSS
        2715 Good BOSS
        2716 Good BOSS
        2717 Good BOSS
        2718 Good BOSS

    """
    def __init__(self, columns, names=None):

        self.session = session

        if not isinstance(columns, list):
            columns = [columns]
        for ii, col in enumerate(columns):
            if isinstance(col, basestring):
                columns[ii] = eval(col)
            if not isinstance(columns[ii], InstrumentedAttribute):
                raise TypeError('{0} is not an InstrumentedAttribute'.format(columns[ii]))

        tables =

        if query is None:
            self.query = self.session.query(Exposure)
        else:
            self.query = query


    def joinTable(self, tables):
        """Joins a table or series of tables to the query."""
        if not isinstance(tables, list):
            tables = [tables]
        for tt in tables:
            self.query = self.query.join(tt)

    def runBatch(self, *args):
        """Allows to run a series of `ExposureQuery` methods in a batch.

        Parameters
        ----------
        args : list of strings
            A series of strings with the names of the `ExposureQuery`
            methods to be called.

        """
        for arg in args:
            if not isinstance(arg, basestring):
                raise TypeError('runBatch attributes must be strings')
            try:
                exec('self.{0}'.format(arg))
            except:
                try:
                    exec('self.{0}()'.format(arg))
                except:
                    raise NameError('{0} is not an PlateDB.queries. '.format(arg) +
                                    'ExposureQuery method.')

    def filterGoodExposures(self):
        """Filters away all bad exposures,"""
        self.query = self.query.filter(Exposure.exposure_status_pk == 1)

    def filterSurvey(self, surveyName):
        """Filters exposures with from a given survey.

        Parameters
        ----------
        surveyName : str
            The name of the survey, to match with platedb.survey.label.

        """
        self.query = self.query.filter(Survey.label.like(surveyName))

    def filterScience(self):
        """Filters only science exposures."""
        self.query = self.query.filter(ExposureFlavor.label.like('Science'))

    def minExposureTime(self, minExpTime=600):
        """Filters away exposures shorter than a certain time.

        Parameters
        ----------
        minExpTime : float
            The minimum exposure time. All rows with exposures shorter
            than this won't be selected.

        """
        self.query = self.query.filter(Exposure.exposure_time >= minExpTime)

    # def saveQuery(self, cols, names=None, toFile=None,
    #               format='ascii.fixed_width', **kwargs):
    #     """Creates an `astropy.table` of the query and saves it to disk.

    #     This method returns an `astropy.table` of the query. If the toFile
    #     options is not None, it will also save it to disk in one of the formats
    #     supported by `astropy.table`.

    #     Parameters
    #     ----------
    #     cols : list of str
    #         A list with the attributes of dataModel.Exposure to be included in the
    #         list.
    #     names : list of str
    #         The column titles for the data.
    #     toFile : None or str
    #         If a string, the name of the file where the data will be saved. If None,
    #         it only returns the table.
    #     format : str
    #         An `astropy.table.write` compatible format.
    #     kwargs
    #         Any other argument that will be pased to `astropy.table.write`.

    #     Examples
    #     --------

    #       Examples ::

    #         qq = ExposureQuery()
    #         tab = qq.saveQuery(cols=['start_time', 'survey.label'],
    #                            names=['StartTime', 'Survey'], toFile='MyFile.dat',
    #                            format='ascii.fixed_width', delimiter=';')

    #     """

    #     import astropy.table as table

    #     if isinstance(cols, list):
    #         cols = ', '.join(cols)

    #     rows = eval('list(self.query.with_entities({0}))'.format(cols))
    #     tab = table.Table(zip(*rows), names=names)

    #     if toFile is not None:
    #         tab.write(toFile, format=format, **kwargs)

    #     return tab
