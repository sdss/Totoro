#!/usr/bin/env python
# encoding: utf-8
"""

footprint.py

Created by José Sánchez-Gallego on 1 Mar 2016.
Licensed under a 3-clause BSD license.

Revision history:
    1 Mar 2016 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib as mpl
from matplotlib import path
from matplotlib.patches import PathPatch
from mangaLib.data import UKIDSS as UKIDSS_Regions


__all__ = ['HSC', 'UKIDSS', 'ATLAS', 'ALFALFA', 'ApertifMedDeep', 'GAMA',
           'CVn', 'HETDEX', 'PerseusPisces', 'getPlatesInFootprint']


def etalambda2radec(eta, lambd):

    racen = 185.0
    deccen = 32.5
    r2d = 180. / np.pi
    d2r = 1. / r2d

    coslambda = np.cos(lambd * d2r)
    dec = r2d * np.arcsin(coslambda * np.sin((eta + deccen) * d2r))
    ra = r2d * np.arctan2(np.sin(lambd * d2r),
                          coslambda * np.cos((eta + deccen) * d2r)) + racen

    return np.hstack((ra[np.newaxis].T, dec[np.newaxis].T))


def getSDSSRegion(vertices):
    """Returns a patch from SDSS coordinates."""

    eta = np.linspace(vertices[0], vertices[1], 5e1)
    lambd = np.linspace(vertices[2], vertices[3], 5e1)

    coords = etalambda2radec(eta, vertices[2])
    coords = np.concatenate((coords, etalambda2radec(vertices[1], lambd)))
    coords = np.concatenate((coords, etalambda2radec(eta[::-1], vertices[3])))
    coords = np.concatenate((coords, etalambda2radec(vertices[0],
                                                     lambd[::-1])))

    regPath = path.Path(coords, closed=True)
    regPatch = PathPatch(regPath)

    return regPatch


def getPolygon(vertices):
    """Returns a patch with polygonal shape based on the input vertices."""

    regPath = path.Path(vertices, closed=True)
    regPatch = PathPatch(regPath)

    return regPatch


def getRectangle(vertices, angle=0.0):
    """Returns a rectangular patch."""

    ra0 = np.min(vertices[0:2])
    ra1 = np.max(vertices[0:2])
    dec0 = np.min(vertices[2:])
    dec1 = np.max(vertices[2:])

    coords = np.array([[ra0, dec0],
                       [ra0, dec1],
                       [ra1, dec1],
                       [ra1, dec0],
                       [ra0, dec0]])
    regPath = path.Path(coords, closed=True)

    recPatch = PathPatch(regPath)

    if angle != 0.0:
        xCentre = 0.5 * (ra0 + ra1)
        yCentre = 0.5 * (dec0 + dec1)
        tt = mpl.transforms.Affine2D().rotate_deg_around(xCentre,
                                                         yCentre,
                                                         angle)

        trans = tt
        recPatch.set_transform(trans)

    return recPatch

# HSC regions
HSC_1 = getRectangle((22 * 15, 360, -1, 7))
HSC_2 = getRectangle((0, 2.6666 * 15, -1, 7))
HSC_3 = getRectangle((1.83333 * 15, 2.66666 * 15, -7, -1))
HSC_4 = getRectangle((8.5 * 15, 15 * 15, -2, 5))
HSC_5 = getRectangle((13.3 * 15, 16.6666 * 15, 42.5, 44))
HSC = [HSC_1, HSC_2, HSC_3, HSC_4, HSC_5]

# A special region which is the HSC-N regions a bit wider
HSC_N_Wide = getRectangle((13.3 * 15, 16.6666 * 15, 41.5, 45))

Stripe82 = getRectangle((0, 360, -1, 1))

SGC = [getRectangle((0, 5 * 15, -20, 20)), getRectangle((300, 360, -20, 20))]

# Herschel-ATLAS region
ATLAS = getRectangle((199.5 - 7.5, 199.5 + 7.5, 29 - 5, 29 + 5), angle=-8)

# GAMA fields
GAMA_1 = getRectangle((129, 141, -1, 3))
GAMA_2 = getRectangle((174, 186, -2, 2))
GAMA_3 = getRectangle((211.5, 223.5, -2, 2))
GAMA = [GAMA_1, GAMA_2, GAMA_3]

# ALFALFA regions
ALFALFA_1 = getRectangle((7.5 * 15, 16.5 * 15, 0, 36))
ALFALFA_2 = getRectangle((0, 3. * 15, 0, 36))
ALFALFA_3 = getRectangle((22 * 15, 24 * 15, 0, 36))
ALFALFA = [ALFALFA_1, ALFALFA_2, ALFALFA_3]

# HETDEX field
HETDEX_vertices = np.array(
    [[164.0, 45.5],
     [180.0, 48.9],
     [196.0, 49.7],
     [210.0, 48.8],
     [226.0, 45.0],
     [229.0, 51.0],
     [210.0, 55.4],
     [196.0, 56.2],
     [180.0, 55.3],
     [160.5, 51.0],
     [164.0, 45.5]])

CVn_vertices = np.array(
    [[184.98, 31.55],
     [191.76, 31.55],
     [192.51, 46.30],
     [184.19, 46.30],
     [184.98, 31.55]])

HETDEX = getPolygon(HETDEX_vertices)
PerseusPisces = getRectangle((23.7, 33.3, 29.9, 37.9))
CVn = getPolygon(CVn_vertices)

# Apertif medium-depth fields
ApertifMedDeep = [HETDEX, PerseusPisces, CVn]

# Latest UKIDSS footprint
UKIDSS = [getPolygon(UKIDSS_Regions.UKIDSS_Region1),
          getPolygon(UKIDSS_Regions.UKIDSS_Region2),
          getPolygon(UKIDSS_Regions.UKIDSS_Region3),
          getPolygon(UKIDSS_Regions.UKIDSS_Region4),
          getPolygon(UKIDSS_Regions.UKIDSS_Region5),
          getPolygon(UKIDSS_Regions.UKIDSS_Region6)]


# Some renaming, for convenience
UKIDSS_240 = UKIDSS[4]
UKIDSS_120 = UKIDSS[5]
UKIDSS_ATLAS = UKIDSS[3]
ALFALFA_NGC = ALFALFA[0]


def getPlatesInFootprint(plates):
    """Returns the list of plates that are within MaNGA footprint."""

    plates = np.atleast_1d(plates)

    footprintPlates = []
    for plate in plates:
        if ((UKIDSS_120.get_path().contains_points([plate.coords])[0] and
            ALFALFA_NGC.get_path().contains_points([plate.coords])[0]) or
           (UKIDSS_240.get_path().contains_points([plate.coords])[0] and
            ALFALFA_NGC.get_path().contains_points([plate.coords])[0]) or
           ATLAS.get_path().contains_points([plate.coords])[0] or
           UKIDSS_ATLAS.get_path().contains_points([plate.coords])[0] or
           SGC[0].get_path().contains_points([plate.coords])[0] or
           SGC[1].get_path().contains_points([plate.coords])[0] or
           HETDEX.get_path().contains_points([plate.coords])[0] or
           PerseusPisces.get_path().contains_points([plate.coords])[0] or
           CVn.get_path().contains_points([plate.coords])[0] or
           HSC_N_Wide.get_path().contains_points([plate.coords])[0]):

            footprintPlates.append(plate)

    if len(footprintPlates) == 1:
        return footprintPlates[0]
    else:
        return footprintPlates
