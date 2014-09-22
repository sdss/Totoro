#!/usr/bin/env python
# encoding: utf-8
"""
timeline.py

Created by José Sánchez-Gallego on 1 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    1 Aug 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.internal.manga.Totoro import config, log, site, TotoroDBConnection
import warnings
from sdss.internal.manga.Totoro.exceptions import TotoroUserWarning
from sdss.internal.manga.Totoro import utils
from sdss.internal.manga.Totoro import logic


db = TotoroDBConnection()
session = db.Session()


class Timelines(list):

    def __init__(self, observingBlocks, plates=None, mode='planner', **kwargs):

        self.mode = mode

        initList = []
        for row in observingBlocks:
            initList.append(Timeline(row['JD0'], row['JD1'], plates=plates))

        list.__init__(self, initList)

    def schedule(self, mode=None, **kwargs):

        if mode is None:
            mode = self.mode

        plates = self[0]._plates

        for timeline in self:
            timeline._plates = plates
            timeline.schedule(mode=mode, **kwargs)
            plates = timeline._plates


class Timeline(object):

    def __init__(self, startTime, endTime, plates=None, **kwargs):

        from sdss.internal.manga.Totoro import dbclasses

        self.startTime = startTime
        self.endTime = endTime

        if plates is None:
            self._plates = dbclasses.Plates()
        else:
            self._plates = plates

        self.site = site

    def schedule(self, mode='planner', **kwargs):

        log.info('scheduling timeline with JD0={0:.4f}, JD1={1:.4f}'
                 .format(self.startTime, self.endTime))

        currentTime = self.startTime
        remainingTime = utils.JDdiff(currentTime, self.endTime)

        expTime = (config['exposure']['exposureTime'] /
                   config[mode]['efficiency'])
        maxLeftoverTime = (config[mode]['maxLeftoverTime'] /
                           config[mode]['efficiency'])

        nCarts = 1
        maxNCarts = config[mode]['nCarts']

        while remainingTime >= maxLeftoverTime and nCarts <= maxNCarts:

            # log.info('Simulating plates for {0:.5f}'.format(currentTime))

            optimalPlate, newExposures = logic.getOptimalPlate(
                self._plates, currentTime, self.endTime, expTime=expTime,
                mode=mode, **kwargs)

            if optimalPlate is None:

                warnings.warn('no valid plates found at JD={0:.4f} '
                              '(timeline ends at JD={1:.4f})'.format(
                                  currentTime, self.endTime),
                              TotoroUserWarning)
                currentTime += config['exposure']['exposureTime'] / 86400.

            else:

                for exp in newExposures:
                    if exp is None:
                        continue
                    optimalPlate.addMockExposure(exposure=exp)

                nCarts += 1
                newTime = optimalPlate.getLastExposure().getJD()[1]
                currentTime = newTime

            remainingTime = utils.JDdiff(currentTime, self.endTime)

        return
