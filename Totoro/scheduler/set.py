#!/usr/bin/env python
# encoding: utf-8
"""
set.py

Created by José Sánchez-Gallego on 18 Mar 2014.
Licensed under a 3-clause BSD license.

Revision history:
    18 Mar 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import numpy as np
from numbers import Real
from astropy.coordinates import Longitude
from astropy import units as uu
from .exposure import Exposure
from ..core.defaults import DITHER_POSITIONS_NEEDED, SN2_FACTOR_SET
from ..core.defaults import SEEING_POOR, SEEING_EXCELLENT, MAX_DIFF_SEEING_SET


class Set(list):

    def __init__(self, ID=None, exposures=None, complete=None, avgSeeing=None,
                 SN_red=None, SN_blue=None, HAlimits=None,
                 visibilityWindow=None, setQuality=None, **kwargs):

        if exposures is not None:
            list.__init__(self, exposures)

        self.ID = ID
        self._complete = complete
        self._avgSeeing = avgSeeing
        self._SN_red = SN_red
        self._SN_blue = SN_blue
        self._HAlimits = HAlimits

        if visibilityWindow is None:
            self.visibilityWindow = visibilityWindow
        elif isinstance(visibilityWindow,
                        (Longitude, tuple, list, np.ndarray)):
            self.visibilityWindow = Longitude(visibilityWindow, unit=uu.hour)
        else:
            raise ValueError('visibilityWindow must be a list.')

        self._setQuality = setQuality
        self._kwargs = kwargs

    def addExposure(self, **kwargs):
        self.append(Exposure(Set=self, **kwargs))

    def getValidExposures(self):

        assert self.visibilityWindow is None or isinstance(
            self.visibilityWindow, (list, tuple, np.ndarray, Longitude))

        if self.visibilityWindow is None:
            newList = [exp for exp in self if exp.valid is True
                       and exp.obsType == 'sci']
        else:
            newList = []
            visWindow = Longitude(self.visibilityWindow, unit=uu.hour)
            for exp in self:
                if (exp.HAstart >= visWindow[0] and
                        exp.HAend <= visWindow[1] and
                        exp.valid is True and
                        exp.obsType == 'sci'):
                    newList.append(exp)

        return Set(ID=self.ID, exposures=newList, complete=self._complete,
                   avgSeeing=self._avgSeeing, SN_red=self._SN_red,
                   SN_blue=self._SN_blue, HAlimits=self._HAlimits,
                   setQuality=self._setQuality, **self._kwargs)

    def _getCompletionStatus(self):

        snOK = False
        ditherPosOK = False
        HAlimitsOK = False
        seeingOK = False

        validExp = self.getValidExposures()

        ditherPos = [exp.ditherPos for exp in validExp]
        if (all([dd in DITHER_POSITIONS_NEEDED() for dd in ditherPos]) and
                all([dd in ditherPos for dd in DITHER_POSITIONS_NEEDED()])):
            ditherPosOK = True

        if (self.get_HA_minmax() is not None and
                Longitude(self.get_HA_minmax()[1] - self.get_HA_minmax()[0],
                          unit=uu.hour)
                <= Longitude('1h')):
            HAlimitsOK = True

        if validExp.nExposures >= 1:
            minMaxSN = self._getMaxMinSN()
            if (all(minMaxSN[1, :]/minMaxSN[0, :] <= SN2_FACTOR_SET()) and
                    all(minMaxSN[3, :]/minMaxSN[2, :] <= SN2_FACTOR_SET())):
                snOK = True

        avgSeeings = [exp.avgSeeing for exp in validExp]
        if np.max(avgSeeings) - np.min(avgSeeings) <= MAX_DIFF_SEEING_SET:
            seeingOK = True

        return (bool(snOK * ditherPosOK * HAlimitsOK),
                (snOK, ditherPosOK, HAlimitsOK, seeingOK))

    def getFailedTests(self):

        tests = self._getCompletionStatus()[1]
        labels = ['SN2', 'ditherPositions', 'HAlimits', 'avgSeeing']
        return [labels[ii] for ii in range(len(labels)) if tests[ii] is False]

    def _getMaxMinSN(self):
        """Returns an array with the max and min SN in the valid exposures.

        The array returned is of the form
            [[minSN_Blue_1, minSN_Blue_2],
             [maxSN_Blue_1, maxSN_Blue_2],
             [minSN_Red_1, minSN_Red_2],
             [maxSN_Red_1, maxSN_Red_2]]

        """

        validExp = self.getValidExposures()

        if validExp.nExposures == 0:
            return np.array([[0., 0.], [0., 0.], [0., 0.], [0., 0.]])

        SN_red = np.array([exp.SN_red for exp in validExp])
        SN_blue = np.array([exp.SN_blue for exp in validExp])
        maxSNred = np.max(SN_red, axis=0)
        minSNred = np.min(SN_red, axis=0)
        maxSNblue = np.max(SN_blue, axis=0)
        minSNblue = np.min(SN_blue, axis=0)

        return np.array([minSNblue, maxSNblue, minSNred, maxSNred])

    def get_HA_minmax(self, exposures=None):

        if exposures is None:
            exposures = self.getValidExposures()

        if len(exposures) == 0:
            return None

        HAs = []
        for exp in exposures:
            HAs += [exp.HAstart.hour, exp.HAend.hour]
        return Longitude([np.min(HAs), np.max(HAs)], unit=uu.hour)

    @property
    def nExposures(self):
        return len(self.getValidExposures())

    @property
    def complete(self):
        if self._complete is not None:
            return self._complete
        return self._getCompletionStatus()[0]

    @complete.setter
    def complete(self, value):
        assert isinstance(value, bool) and value is not None
        self._complete = value

    @property
    def avgSeeing(self):
        if self._avgSeeing is not None:
            return self._avgSeeing
        return np.mean([exp.avgSeeing for exp in self.getValidExposures()])

    @avgSeeing.setter
    def avgSeeing(self, value):
        assert isinstance(value, Real)
        self._avgSeeing = value

    @property
    def SN_red(self):
        if self._SN_red is not None:
            return self._SN_red
        if self.nExposures == 0:
            return np.array([0.0, 0.0])
        return np.sum(np.atleast_2d(
            [exp.SN_red for exp in self.getValidExposures()]), axis=0)

    @SN_red.setter
    def SN_red(self, value):
        assert isinstance(value, (np.ndarray, list, tuple))
        self._SN_red = value

    @property
    def SN_blue(self):
        if self._SN_blue is not None:
            return self._SN_blue
        if self.nExposures == 0:
            return np.array([0.0, 0.0])
        return np.sum(np.atleast_2d(
            [exp.SN_blue for exp in self.getValidExposures()]), axis=0)

    @SN_blue.setter
    def SN_blue(self, value):
        assert isinstance(value, (np.ndarray, list, tuple))
        self._SN_blue = value

    @property
    def HAlimits(self):
        if self._HAlimits is not None:
            return self._HAlimits

        HAmaxmin = self.get_HA_minmax()
        if HAmaxmin is None:
            if self.visibilityWindow is not None:
                return self.visibilityWindow
            return None

        oneHour = Longitude('1h')
        if Longitude(HAmaxmin[1] - HAmaxmin[0], unit=uu.hour) > oneHour:
            return HAmaxmin

        return Longitude(
            [Longitude(HAmaxmin[1]-oneHour), Longitude(HAmaxmin[0]+oneHour)])

    @HAlimits.setter
    def HAlimits(self, value):
        assert isinstance(value, (np.ndarray, list, tuple))
        self._HAlimits = np.array(value)

    @property
    def setQuality(self):
        if self._setQuality is not None:
            return self._setQuality

        if self.complete is False:
            return 'bad'

        if self.avgSeeing > SEEING_POOR():
            return 'poor'
        elif self.avgSeeing <= SEEING_POOR() and \
                self.avgSeeing > SEEING_EXCELLENT():
            return 'good'
        else:
            return 'excellent'

    @setQuality.setter
    def setQuality(self, value):
        assert isinstance(value, basestring)
        value = value.lower()
        self._setQuality = value
        if value in ['good', 'excellent', 'poor']:
            self.complete = True
        elif value == 'bad':
            self.complete = False
        else:
            raise ValueError('setQuality cannot have value {0}'.format(value))
