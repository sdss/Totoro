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

UKIDSS_Region1 = np.array([
    [061.95838, +00.87955],
    [054.13354, +06.28934],
    [48.88450, +06.95682],
    [046.32311, +08.43322],
    [039.96898, +08.81795],
    [035.20637, +13.01168],
    [035.12406, +14.31154],
    [030.93135, +17.09270],
    [0.0000000, +17.00000],
    [0.0000000, -02.00000],
    [027.91088, -02.19216],
    [047.03013, -02.55489],
    [059.19465, -02.59425],
    [062.08225, -00.05170],
    [061.95838, +00.87955]])

UKIDSS_Region2 = np.array([
    [360.00000, +17.00000],
    [351.87501, +16.96904],
    [349.93036, +17.04856],
    [345.16645, +15.03193],
    [339.41813, +11.29414],
    [332.81999, +05.71610],
    [326.77426, +02.16565],
    [308.18905, +02.31217],
    [308.04936, -02.10125],
    [333.13159, -02.30051],
    [359.89682, -02.32808],
    [360.00000, +17.00000]])

UKIDSS_Region3 = np.array([
    [239.36011, -02.72270],
    [239.49965, +02.40468],
    [240.41374, +09.81649],
    [237.14271, +11.35246],
    [220.12234, +15.07941],
    [209.19142, +15.81317],
    [206.64129, +17.20587],
    [180.21959, +16.95451],
    [178.98269, +17.89222],
    [167.77638, +17.63917],
    [166.22939, +16.14737],
    [155.82172, +15.56455],
    [148.06986, +14.87760],
    [123.50277, +09.32339],
    [124.93499, -02.63832],
    [133.96699, -04.23402],
    [169.63894, -03.19124],
    [171.64508, -05.06182],
    [186.01105, -04.97972],
    [205.37750, -04.99526],
    [207.43854, -03.45983],
    [225.76219, -03.65258],
    [232.18386, -03.77643],
    [239.36011, -02.72270]])

UKIDSS_Region4 = np.array([
    [211.18558, +21.91363],
    [211.38403, +37.05334],
    [188.58584, +36.73161],
    [189.76090, +21.58410],
    [198.23781, +20.75589],
    [211.18558, +21.91363]])

UKIDSS_Region5 = np.array([
    [242.50593, +20.93451],
    [230.69122, +27.35972],
    [220.47165, +32.45357],
    [222.59255, +35.00932],
    [226.67155, +34.71249],
    [232.46691, +32.07994],
    [234.40908, +31.63117],
    [238.06585, +28.74299],
    [240.20232, +31.69502],
    [244.57836, +32.74295],
    [252.16101, +33.78028],
    [255.60497, +30.72946],
    [250.29212, +25.29416],
    [245.53296, +21.25760],
    [242.50593, +20.93451]])

UKIDSS_Region6 = np.array([
    [120.38225, +17.00487],
    [114.65075, +20.61480],
    [110.51386, +25.41934],
    [112.14745, +27.20383],
    [111.99926, +28.35569],
    [114.62997, +30.68592],
    [132.00404, +30.98818],
    [135.50603, +33.43509],
    [136.85960, +34.49213],
    [143.02820, +36.78138],
    [147.04215, +36.75095],
    [149.45191, +34.14314],
    [147.78351, +31.55026],
    [145.65717, +29.22833],
    [138.61797, +26.43600],
    [130.26806, +22.78999],
    [125.86197, +19.94148],
    [121.75807, +17.33957],
    [120.38225, +17.00487]])

UKIDSS = [getPolygon(UKIDSS_Region1),
          getPolygon(UKIDSS_Region2),
          getPolygon(UKIDSS_Region3),
          getPolygon(UKIDSS_Region4),
          getPolygon(UKIDSS_Region5),
          getPolygon(UKIDSS_Region6)]


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
