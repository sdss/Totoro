#!/usr/bin/env python
# encoding: utf-8
"""
exposure.py

Created by José Sánchez-Gallego on 13 Mar 2014.
Licensed under a 3-clause BSD license.

Revision history:
    13 Mar 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from ..exceptions import TotoroError, NoMangaExposure
from .. import TotoroDBConnection
from ..logic import mangaLogic
from .. import log, config, site
import numpy as np
from .. import utils
from astropy import time
from .. import dustMap
import warnings


totoroDB = TotoroDBConnection()
plateDB = totoroDB.plateDB
mangaDB = totoroDB.mangaDB
session = totoroDB.Session()


class Exposure(plateDB.Exposure):

    def __new__(cls, input=None, format='pk', parent='plateDB', *args,
                **kwargs):

        if input is None:
            return plateDB.Exposure.__new__(cls)

        base = cls.__bases__[0]

        if isinstance(input, base):
            instance = input

        elif isinstance(input, mangaDB.Exposure):
            instance = input.platedbExposure

        else:
            if parent.lower() == 'platedb':
                with session.begin(subtransactions=True):
                    instance = session.query(base).filter(
                        eval('{0}.{1} == {2}'.format(base.__name__,
                                                     format, input))).one()

            elif parent.lower() == 'mangadb':
                with session.begin(subtransactions=True):
                    instance = session.query(base).join(
                        mangaDB.Exposure).filter(
                            eval('mangaDB.Exposure.{0} == {1}'.format(
                                format, input))).one()

        instance.__class__ = cls

        return instance

    def __init__(self, input, format='pk', mock=False, silent=False,
                 *args, **kwargs):

        self._valid = None
        self._ditherPosition = None
        self._sn2Array = None
        self._seeing = None

        self.isMock = mock
        self.kwargs = kwargs

        if not silent:
            log.debug('Loaded exposure plateDB.Exposure.pk={0}'
                      .format(self.pk))

        self.mlhalimit = utils.mlhalimit(self.dec)

        self._mangaExposure = (self.mangadbExposure[0]
                               if len(self.mangadbExposure) > 0 else
                               mangaDB.Exposure())
        if self._mangaExposure is None:
            warnings.warn('plateDB.Exposure.pk={0} has no mangaDB.Exposure '
                          'counterpart.', NoMangaExposure)

    def __repr__(self):
        return ('<Totoro Exposure (mangaDB.Exposure.pk={0}, exposure_no={1}, '
                'ditherPos={2})>'
                .format(self._mangaExposure.pk, self.exposure_no,
                        self.ditherPosition))

    def update(self, **kwargs):
        """Reloads the exposure."""

        newSelf = Exposure(self.pk, fromat='pk', silent=True, mock=self.isMock,
                           **self.kwargs)
        self = newSelf

        log.debug('Exposure pk={0} has been reloaded'.format(
                  self.pk))

    @classmethod
    def createMockExposure(cls, startTime=None, expTime=None,
                           ditherPosition=None, ra=None, dec=None,
                           silent=False, **kwargs):

        if ra is None or dec is None:
            raise TotoroError('ra and dec must be specified')

        newExposure = Exposure(None, mock=True, ra=ra, dec=dec, **kwargs)
        newExposure.pk = None if 'pk' not in kwargs else kwargs['pk']

        if startTime is None:
            startTime = time.Time.now().jd

        if expTime is None:
            expTime = config['exposure']['exposureTime']

        newExposure.ditherPosition = ditherPosition

        tt = time.Time(startTime, format='jd', scale='tai')
        t0 = time.Time(0, format='mjd', scale='tai')
        startTimePlateDB = (tt - t0).sec
        newExposure.start_time = startTimePlateDB

        newExposure.exposure_time = expTime

        newExposure.simulateObservedParamters()

        if not silent:
            log.debug('Created mock exposure with ra={0:.3f} and dec={0:.3f}'
                      .format(newExposure.ra, newExposure.dec))

        return newExposure

    def simulateObservedParamters(self):

        self._seeing = 1.0

        self._dust = dustMap(self.ra, self.dec)

        haRange = self.getHA()
        ha = utils.calculateMean(haRange)
        self._airmass = utils.computeAirmass(self.dec, ha)

        sn2Red = config['simulation']['redSN2'] / self._airmass ** \
            config['simulation']['alphaRed'] / self._dust['iIncrease'][0]
        sn2Blue = config['simulation']['blueSN2'] / self._airmass ** \
            config['simulation']['alphaBlue'] / self._dust['gIncrease'][0]

        self._sn2Array = np.array([sn2Blue, sn2Blue, sn2Red, sn2Red])

        # self._valid = True
        # self.status = 'Good'

    @property
    def ra(self):
        return self.getCoordinates()[0]

    @property
    def dec(self):
        return self.getCoordinates()[1]

    def getCoordinates(self):

        if 'ra' in self.kwargs and 'dec' in self.kwargs:
            if (self.kwargs['ra'] is not None and
                    self.kwargs['dec'] is not None):
                return np.array(
                    [self.kwargs['ra'], self.kwargs['dec']], np.float)
        else:
            pointing = (self.observation.plate_pointing.pointing)
            return np.array(
                [pointing.center_ra, pointing.center_dec], np.float)

    def getHA(self):
        """Returns the HA range in which the exposure was taken."""

        startTime = float(self.start_time)
        expTime = float(self.exposure_time)

        t0 = time.Time(0, format='mjd', scale='tai')
        tStart = t0 + time.TimeDelta(startTime, format='sec', scale='tai')

        lst = site.localSiderialTime(tStart.jd)
        ha = (lst * 15. - self.ra) % 360.

        return np.array([ha, ha + expTime / 3600. * 15]) % 360.

    def getSN2Array(self):
        """Returns an array with the SN2 of the exposure. The return
        format is [b1SN2, b2SN2, r1SN2, r2SN2]."""

        if self._sn2Array is not None:
            return self._sn2Array

        SN2Values = self._mangaExposure.sn2values[0]

        return np.array([SN2Values.b1_sn2, SN2Values.b2_sn2,
                         SN2Values.r1_sn2, SN2Values.r2_sn2])

    @property
    def valid(self):
        if self._valid is not None:
            return self._valid
        else:
            return self.isValid()[0]

    @valid.setter
    def valid(self, value):
        self._valid = value

    def isValid(self, flag=True, **kwargs):
        """Checks if an exposure is valid."""

        return mangaLogic.checkExposure(self, flag=flag, **kwargs)

    @property
    def ditherPosition(self):
        if self._ditherPosition is None and self.isMock is False:
            return self._mangaExposure.dither_position[0].upper()
        else:
            return self._ditherPosition

    @ditherPosition.setter
    def ditherPosition(self, value):
        self._ditherPosition = value

    @property
    def seeing(self):
        if not self.isMock and self._seeing is None:
            return self._mangaExposure.seeing
        else:
            return self._seeing

    def getLST(self):

        ha0, ha1 = self.getHA()

        lst0 = (ha0 + self.ra) % 360. / 15
        lst1 = (ha1 + self.ra) % 360. / 15

        return np.array([lst0, lst1])

    def getUT(self, format=None):

        startTime = float(self.start_time)
        t0 = time.Time(0, format='mjd', scale='tai')

        tStart = t0 + time.TimeDelta(startTime, format='sec', scale='tai')
        tEnd = tStart + time.TimeDelta(float(self.exposure_time), format='sec',
                                       scale='tai')

        ut0 = tStart.datetime
        ut1 = tEnd.datetime

        if format == 'str':
            return ('{0:%H:%M}'.format(ut0), '{0:%H:%M}'.format(ut1))
        else:
            return (ut0, ut1)

        return (ut0, ut1)

    def getJD(self):

        startTime = float(self.start_time)
        t0 = time.Time(0, format='mjd', scale='tai')

        tStart = t0 + time.TimeDelta(startTime, format='sec', scale='tai')
        tEnd = tStart + time.TimeDelta(
            float(self.exposure_time), format='sec', scale='tai')

        return (tStart.jd, tEnd.jd)

    def getPlatePK(self):
        return int(self.observation.plate_pointing.plate.pk)

    def getMJD(self):
        return self.mjd()
