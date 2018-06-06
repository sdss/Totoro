#!/usr/bin/env python
# encoding: utf-8
"""
dustMap.py

Created by José Sánchez-Gallego on 13 May 2014.
Licensed under a 3-clause BSD license.

Revision history:
    13 May 2014 J. Sánchez-Gallego
      Initial version

"""

import glob
import os

import numpy as np
from astropy.wcs import WCS
from scipy.interpolate import griddata

from .ccmUnred import ccmUnred


try:
    import fitsio
    fitsLib = 'fitsio'
except:
    try:
        from astropy.io import fits
        fitsLib = 'fits'
    except:
        raise ImportError('either fitsio or astropy.io.fits are needed.')

try:
    from astropysics import coords
    astropysicsLib = True
except:
    try:
        from astropy.coordinates import SkyCoord
        astropysicsLib = False
    except:
        raise ImportError('either astropysics or astropy are needed.')

try:
    __defaultMaps__ = glob.glob(
        os.path.join(os.environ['DUST_DIR'], 'maps', 'SFD_dust_4096_*.fits'))
except:
    __defaultMaps__ = None


class DustMap(object):
    """A class to analyse dust maps and derive values.

    This class allows to load dust maps and determine extinction values at
    specific coordinates. Any dust map with units E(B-V) can be used as long
    as valid WCS information is stored in the header of the FITS file. Both FK5
    or Galactic coordinates can be used.

    Once the class is instantialised, it can be called by providing a list of
    coordinates. The returned dictionary includes the E(B-V) value(s) for the
    input coordinate(s), as well as the iIncrease and gIncrease values at the
    position. iIncrease is defined as the reddening correction for an unitary
    flux in the SDSS i-band:

    iIncrease = 1 / iOut ** 2

    where iOut is the reddened unitary flux in the i-band at the position.

    The returned values are masked arrays with elements masked if the requested
    coordinate was outside the footprint of all of the dust maps.

    Parameters
    ----------
    maps : string, list of strings, optional
        The path or list of paths for the dust maps to be used. If not
        specified, dust maps are searched at the path $DUST_DIR/maps/.

    Example
    -------
    The class can be initialised as ::
      >> from sdss.manga import DustMap
      >> dMap = DustMap()

    And then called as ::
      >> dustValues = dMap(12.5, 30)

    This will return the dust parameters at the position RA=12.5 deg,
    Dec=30 deg. It's possible to call the DustMap instance with a list of
    coordinates. To determine the dust parameters at (12.5, 30) and
    (210, -15) ::
      >> dustValues = dMap([12.5, 210], [30, -15])

    By default, equatorial coordinates are assumed. Galactic coordinates can be
    used by adding the keyword `coordSys='galactic'`. The normal behaviour is
    not to interpolate the values of the dust map. This can be overridden by
    defining `interpolate=True`. The output dictionary is of the form ::
      >> print(dustValues)
      {'EBV': masked_array(data = [0.05857369377550144 0.08294817811791555],
              mask = [False False],
        fill_value = 1e+20),
 'gIncrease': masked_array(data = [1.5123812061910469 1.7964882075835737],
              mask = [False False],
        fill_value = 1e+20),
 'galCoords': [[122.56131852693547, -32.870627753783275],
  [326.54366100057314, 44.705445415210264]],
 'iIncrease': masked_array(data = [1.257851242757 1.383847964727236],
              mask = [False False],
        fill_value = 1e+20)}

    """

    def __init__(self, maps=__defaultMaps__, **kwargs):

        if maps is None:
            raise ValueError('no maps defined and $DUST_DIR is not defined.')

        self.maps = [maps] if not hasattr(maps, '__getitem__') else maps

        if False in map(os.path.exists, self.maps):
            raise NameError('dust map(s) does not exist.')

        if len(self.maps) == 0:
            raise ValueError('no dust maps available.')

        self.loadMaps()

    def __call__(self, xx, yy, coordSys='fk5', interpolate=False):

        if coordSys not in ['fk5', 'galactic']:
            raise ValueError('coordSys must be fk5 or galactic.')

        xx = np.atleast_1d(xx)
        yy = np.atleast_1d(yy)
        if len(xx) != len(yy):
            raise ValueError('x and y coordinates must be of the same length.')

        if coordSys == 'galactic':
            coordsGal = np.array([xx, yy]).T
        else:
            if astropysicsLib:
                ccICRS = [coords.ICRSCoordinates(xx[ii], yy[ii]) for ii in range(len(xx))]
                ccGal = [cc.convert(coords.GalacticCoordinates) for cc in ccICRS]
                coordsGal = np.array([[cc.l.degrees for cc in ccGal],
                                      [cc.b.degrees for cc in ccGal]]).T
            else:
                ccSkyCoords = SkyCoord(xx, yy, frame='icrs', unit='deg')
                ccGal = ccSkyCoords.galactic
                coordsGal = np.array([ccGal.l.deg, ccGal.b.deg]).T

        return self._getValues(coordsGal, interpolate=interpolate)

    def _getValues(self, cc, interpolate=True):

        # Central wavelengths for u, g, r, i, z
        lambdaIn = np.array([3551., 4686., 6165., 7481., 8931.])
        # Unitary flux
        fluxIn = np.array([1., 1., 1., 1., 1.])

        ebv = np.ma.zeros(cc.shape[0])
        ebv.mask = True
        iIncrease = ebv.copy()
        gIncrease = ebv.copy()

        for ii in range(len(self._hdu)):

            ww = self._wcs[ii]
            data = self._data[ii]

            pixels = ww.wcs_world2pix(cc, 0)

            for jj in range(len(pixels)):

                if not np.ma.is_masked(ebv.mask[jj]) and ebv[jj] != 0.0:
                    continue

                xx, yy = pixels[jj]
                x0 = int(xx)
                y0 = int(yy)

                # If the pixel coordinates are within the footprint of the
                # image, stores the value of E(B-V). This value can be
                # overwritten later if interpolate=True, but this ensures that
                # there is a stored value if interpolate=True but the boundary
                # conditions are invalid (the pixel is too close to the border
                # of the image).
                if (x0 >= 0.0 and y0 >= 0.0 and x0 <= data.shape[1] - 1 and
                        y0 <= data.shape[0] - 1):
                    ebv[jj] = data[y0, x0]

                if interpolate is True:
                    # Checks that the boundary conditions are ok for the
                    # interpolation.
                    if (xx < 1.0 or yy < 1.0 or xx > data.shape[1] - 1 or yy > data.shape[0] - 1):
                        continue

                    mgrid = np.mgrid[y0 - 1:y0 + 2, x0 - 1:x0 + 2]
                    grid = np.array([mgrid[0].flatten(), mgrid[1].flatten()]).T
                    slice = data[y0 - 1:y0 + 2, x0 - 1:x0 + 2].flatten()
                    ebv[jj] = griddata(grid, slice, (yy, xx))

            # del data

        # Calculates iIncrease and gIncrease from E(B-V)
        for ii in range(len(ebv)):

            if np.ma.is_masked(ebv[ii]) is True:
                continue

            # Calculates the reddened flux (we use -ebv[ii]).
            fluxOut = ccmUnred(lambdaIn, fluxIn, -ebv[ii])
            iIncrease[ii] = 1. / fluxOut[3]**2
            gIncrease[ii] = 1. / fluxOut[1]**2

        return {
            'galCoords': cc.tolist(),
            'EBV': ebv,
            'iIncrease': iIncrease,
            'gIncrease': gIncrease
        }

    def loadMaps(self):

        self._hdu = []
        self._header = []
        self._wcs = []
        self._data = []

        for mm in self.maps:

            if fitsLib == 'fitsio':
                hduList = fitsio.FITS(mm)
                self._hdu.append(hduList[-1])

                header = hduList[-1].read_header()
                self._header.append(header)
                self._wcs.append(WCS(header))
                self._data.append(hduList[-1].read())

            else:
                hduList = fits.open(mm)
                self._hdu.append(hduList[-1])

                header = hduList[-1].header
                self._header.append(header)
                self._wcs.append(WCS(header))
                self._data.append(hduList[-1].data)
