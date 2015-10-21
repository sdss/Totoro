#!/usr/bin/env python
# encoding: utf-8
"""
plugging.py

Created by José Sánchez-Gallego on 24 Mar 2014.
Licensed under a 3-clause BSD license.

Revision history:
    24 Mar 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from .set import Set
from ..core.defaults import B_SN2, R_SN2
import numpy as np
from ..core.defaults import *
from .. import log


class Plugging(list):

    _complete = False

    def __init__(self, parent=None, sets=[], **kwargs):

        self.parent = parent

        list.__init__(self, sets)

    def addSet(self, **kwargs):

        ids = [ss.ID for ss in self]
        if len(ids) == 0:
            ID = 1
        else:
            ID = np.max(ids) + 1

        log.info('Creating new set ID={0}.'.format(ID))
        ss = Set(ID=ID, plugging=self)
        endTime = ss.addExposures(**kwargs)
        self.append(ss)

        return endTime

    @property
    def complete(self):
        if self._complete is not None:
            return self._complete

        goodSets = [ss for ss in self
                    if ss.quality in ['good', 'excellent']]

        SN_blue = np.sum([ss.SN_blue for ss in goodSets], axis=0)
        SN_red = np.sum([ss.SN_red for ss in goodSets], axis=0)

        if all(SN_red >= R_SN2) and all(SN_blue >= B_SN2):
            return True
        else:
            return False

    @complete.setter
    def complete(self, value):
        assert isinstance(value, bool) or value is None
        self._complete = value

    def hasIncompleteSet(self):
        if any([ss.missingDither for ss in self]):
            return True
        return False

    def getIncompleteSets(self):
        return [ss for ss in self if ss.missingDither]

    def getHAlimit(self):
        if not self.hasIncompleteSet():
            return None
        else:
            return [ss.HAlimits for ss in self if ss.complete is False]

    @property
    def SN_Red(self):
        snRed = np.array([0., 0.])
        for ss in self:
            if ss.quality == 'bad':
                continue
            for ee in ss:
                if ee.valid is True:
                    snRed += ee.SN_red
        return snRed

    @property
    def SN_Blue(self):
        snBlue = np.array([0., 0.])
        for ss in self:
            if ss.quality == 'bad':
                continue
            for ee in ss:
                if ee.valid is True:
                    snBlue += ee.SN_blue
        return snBlue
