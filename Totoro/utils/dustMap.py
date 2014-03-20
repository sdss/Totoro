#!/usr/bin/env python
# encoding: utf-8
"""
dustMapUtils.py

Created by José Sánchez-Gallego on 11 Mar 2014.
Licensed under a 3-clause BSD license.

Revision history:
    11 Mar 2014 J. Sánchez-Gallego
      Initial version

"""


from .ccmUnred import ccmUnred
from astropy.io import fits
import numpy as np
from astropy import coordinates as coo
from astropy import units as uu
from astropy import wcs
import os
from ..exceptions import TotoroError
from ..core.defaults import INCREASE_MAPS, INCREASE_MAPS_FORMAT
from scipy.interpolate import RectBivariateSpline
from numbers import Real


class DustMap(object):

    def __init__(self, maps=INCREASE_MAPS(), format=INCREASE_MAPS_FORMAT(),
                 **kwargs):

        self._files = maps

        if format == 'map':
            if not hasattr(self._files, '__getitem__'):
                self._files = [self._files]
            self._fromMaps(**kwargs)
        elif format == 'grid':
            self._fromGrid(**kwargs)
        else:
            raise TotoroError('format nor understood.')

        self.initGrid()

    def eval(self, xx, yy):
        if isinstance(xx, Real) and isinstance(yy, Real):
            return (self.gIncreaseSpline.ev(xx, yy)[0],
                    self.iIncreaseSpline.ev(xx, yy)[0])
        elif hasattr(xx, '__getitem__') and \
                hasattr(yy, '__getitem__') and len(xx) == len(yy):
            return (self.gIncreaseSpline.ev(xx, yy),
                    self.iIncreaseSpline.ev(xx, yy))
        else:
            raise TotoroError(
                'input must be scalar or arrays of the same size. ')

    def initGrid(self):

        xSize, ySize = self.iIncrease.shape
        xx = np.linspace(self.ra0, self.ra1, xSize)
        yy = np.linspace(self.dec0, self.dec1, ySize)

        self.iIncreaseSpline = RectBivariateSpline(
            xx, yy, self.iIncrease, kx=1, ky=1)

        self.gIncreaseSpline = RectBivariateSpline(
            xx, yy, self.gIncrease, kx=1, ky=1)

    def _fromGrid(self, **kwargs):

        hdu = fits.open(self._files)

        self.gIncrease = hdu[1].data
        self.iIncrease = hdu[2].data

        self.ra0 = hdu[0].header['RA_0']
        self.ra1 = hdu[0].header['RA_1']
        self.dec0 = hdu[0].header['DEC_0']
        self.dec1 = hdu[0].header['DEC_1']

    def _fromMaps(self, ra=[0.0, 360.], dec=[-30, 80.], step=1, **kwargs):

        lambdaIn = np.array([3551., 4686., 6165., 7481., 8931.])
        fluxIn = np.array([1., 1., 1., 1., 1.])

        maps = self._files

        gridRA, gridDec = np.mgrid[ra[0]:ra[1]+step:step,
                                   dec[0]:dec[1]+step:step]

        iIncrease = np.zeros(gridRA.shape, np.float).flatten()
        gIncrease = np.zeros(gridRA.shape, np.float).flatten()
        ebv = np.zeros(gridRA.shape, np.float).flatten()

        icrs = coo.ICRS(gridRA.flatten(), gridDec.flatten(),
                        unit=(uu.degree, uu.degree))

        for mm in maps:

            hdu = fits.open(mm)
            header = hdu[0].header
            data = hdu[0].data
            shape = data.shape
            ww = wcs.WCS(header)

            coords = np.array(
                [icrs.galactic.l.deg, icrs.galactic.b.deg]).T
            pix = ww.wcs_world2pix(coords, 0)
            pix[:, [0, 1]] = pix[:, [1, 0]]
            pix = np.array(pix, dtype=int)

            outOfFramePix = np.where(
                (pix[:, 0] - 1 < 0.) | (pix[:, 0] + 1 >= shape[0]) |
                (pix[:, 1] - 1 < 0.) | (pix[:, 1] + 1 >= shape[1])
            )

            pix[outOfFramePix[0]] = (0.0, 0.0)
            tmpEBV = data[(pix[:, 0], pix[:, 1])]
            tmpEBV[outOfFramePix[0]] = 0.0

            ebv[ebv == 0.0] = tmpEBV[ebv == 0.0]

        for ii in range(len(ebv)):
            fluxOut = ccmUnred(lambdaIn, fluxIn, -ebv[ii])
            iIncrease[ii] = 1. / fluxOut[3] ** 2
            gIncrease[ii] = 1. / fluxOut[1] ** 2

        iIncrease = iIncrease.reshape(gridRA.shape)
        gIncrease = gIncrease.reshape(gridRA.shape)

        self.iIncrease = iIncrease
        self.gIncrease = gIncrease
        self.ra0 = ra[0]
        self.ra1 = ra[1]
        self.dec0 = dec[0]
        self.dec1 = dec[1]

    def save(self, file):
        """Saves the grid as a FITS file."""

        hduPrimary = fits.PrimaryHDU()
        hduPrimary.header.update(
            [('RA_0', self.ra0),
             ('RA_1', self.ra1),
             ('Dec_0', self.dec0),
             ('Dec_1', self.dec1)])

        hduIIncrease = fits.ImageHDU(data=self.iIncrease)
        hduIIncrease.header['EXTNAME'] = 'i_Increase'
        hduGIncrease = fits.ImageHDU(data=self.gIncrease)
        hduGIncrease.header['EXTNAME'] = 'g_Increase'

        hduList = fits.HDUList(
            [hduPrimary, hduGIncrease, hduIIncrease])

        if os.path.exists(file):
            os.remove(file)

        hduList.writeto(file)
