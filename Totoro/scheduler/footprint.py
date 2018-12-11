#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2016-05-01
# @Filename: footprint.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# @Copyright: José Sánchez-Gallego

from __future__ import division, print_function

import matplotlib as mpl
import numpy as np
from matplotlib import path
from matplotlib import pyplot as plt
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Ellipse, PathPatch


__all__ = [
    'HSC', 'UKIDSS', 'ATLAS', 'ALFALFA', 'ApertifMedDeep', 'GAMA', 'CVn', 'HETDEX',
    'PerseusPisces', 'getPlatesInFootprint', 'ALFALFA_Detailed'
]

defaultColours = {
    'GAMA': 'DarkOrchid',
    'HSC': 'k',
    'UKIDSS': 'SaddleBrown',
    'Apertif': 'OliveDrab',
    'ALFALFA': 'b',
    'ALFALFA_Detailed': 'b',
    'ATLAS': 'MediumVioletRed'
}

defaultLW = 1.5


def getAxes(projection='rect'):
    """Returns axes for a particular projection."""

    if projection == 'rect':
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)

        ax.set_xlabel(r'$\alpha_{2000}$')
        ax.set_ylabel(r'$\delta_{2000}$')

        ax.set_xlim(360, 0)
        ax.set_ylim(-20, 80)

    elif projection == 'mollweide':
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111, projection='mollweide')
        org = 120

        tick_labels = np.array([150., 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
        tick_labels = np.remainder(tick_labels + 360 + org, 360)
        tick_labels = np.array(tick_labels / 15., int)

        tickStr = []
        for tick_label in tick_labels[1::2]:
            tickStr.append('')
            tickStr.append('${0:d}^h$'.format(tick_label))

        ax.set_xticklabels(tickStr)  # we add the scale on the x axis
        ax.grid(True)

        ax.set_xlabel(r'$\alpha_{2000}$')
        ax.set_ylabel(r'$\delta_{2000}$')

    return fig, ax


def addLegend(ax, handles, labels, **kwargs):
    """Add a legend using ellipses for the scatter plots."""

    ax.legend(
        handles=handles,
        labels=labels,
        handler_map={Ellipse: HandlerPatch(patch_func=make_legend_ellipse)},
        **kwargs)

    return ax


def make_legend_ellipse(legend, orig_handle, xdescent, ydescent, width, height, fontsize):

    pp = Ellipse(
        xy=(0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent),
        width=(height + ydescent),
        height=(height + ydescent),
        edgecolor='None',
        lw=0.0)

    return pp


def plotEllipse(ax, RA, Dec, org=None, size=3.0, bgcolor='b', zorder=0, alpha=0.8):

    if org:
        RA = np.remainder(RA + 360 - org, 360)  # shift RA values
        ind = RA > 180.
        RA[ind] -= 360  # scale conversion to [-180, 180]
        RA = -RA  # reverse the scale: East to the left

    for ii in range(len(RA)):
        ell = Ellipse(
            xy=(RA[ii], Dec[ii]),
            width=size / np.cos(np.radians(Dec[ii])),
            height=size,
            edgecolor='None',
            facecolor=bgcolor,
            zorder=zorder,
            alpha=alpha,
            lw=0.0)
        ax.add_patch(ell)

    return ell


def etalambda2radec(eta, lambd):

    racen = 185.0
    deccen = 32.5
    r2d = 180. / np.pi
    d2r = 1. / r2d

    coslambda = np.cos(lambd * d2r)
    dec = r2d * np.arcsin(coslambda * np.sin((eta + deccen) * d2r))
    ra = r2d * np.arctan2(np.sin(lambd * d2r), coslambda * np.cos((eta + deccen) * d2r)) + racen

    return np.hstack((ra[np.newaxis].T, dec[np.newaxis].T))


def getSDSSRegion(vertices):
    """Returns a patch from SDSS coordinates."""

    eta = np.linspace(vertices[0], vertices[1], 5e1)
    lambd = np.linspace(vertices[2], vertices[3], 5e1)

    coords = etalambda2radec(eta, vertices[2])
    coords = np.concatenate((coords, etalambda2radec(vertices[1], lambd)))
    coords = np.concatenate((coords, etalambda2radec(eta[::-1], vertices[3])))
    coords = np.concatenate((coords, etalambda2radec(vertices[0], lambd[::-1])))

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

    coords = np.array([[ra0, dec0], [ra0, dec1], [ra1, dec1], [ra1, dec0], [ra0, dec0]])
    regPath = path.Path(coords, closed=True)

    recPatch = PathPatch(regPath)

    if angle != 0.0:
        xCentre = 0.5 * (ra0 + ra1)
        yCentre = 0.5 * (dec0 + dec1)
        tt = mpl.transforms.Affine2D().rotate_deg_around(xCentre, yCentre, angle)

        trans = tt
        recPatch.set_transform(trans)

    return recPatch


def plotPatch(ax, regPatch, zorder=100, projection='rect', useRadians=False, org=0, **kwargs):

    vertices = regPatch.get_path().vertices.copy()

    if projection != 'rect':
        RA = vertices[:, 0]
        RA = np.remainder(RA + 360 - org, 360)  # shift RA values
        ind = RA > 180.
        RA[ind] -= 360  # scale conversion to [-180, 180]
        RA = -RA  # reverse the scale: East to the left

        vertices[:, 0] = RA

    if useRadians:
        vertices = vertices * np.pi / 180.

    transformedPath = path.Path(vertices)

    lw = kwargs.get('lw', defaultLW)

    pathPatch = PathPatch(
        transformedPath,
        edgecolor=kwargs['color'],
        lw=lw,
        facecolor='None',
        alpha=1.,
        zorder=zorder)

    ax.add_patch(pathPatch)

    return


def addText(ax,
            xx,
            yy,
            text,
            ha='left',
            size=None,
            projection='rect',
            useRadians=False,
            org=0,
            **kwargs):

    if projection != 'rect':

        if not size:
            size = 10

        xx = np.remainder(xx + 360 - org, 360)  # shift RA values
        if xx > 180.:
            xx -= 360  # scale conversion to [-180, 180]
        xx = -xx  # reverse the scale: East to the left

    else:

        if not size:
            size = 12

    if useRadians:
        xx = xx * np.pi / 180.
        yy = yy * np.pi / 180.

    ax.text(
        xx,
        yy,
        text,
        size=size,
        color=kwargs['color'],
        fontdict={'family': 'sans-serif'},
        alpha=1.0,
        zorder=200,
        weight='heavy',
        ha=ha)
    mpl.rc(mpl.rcParamsOrig)

    return


def plotRectangle(ax,
                  regPatch,
                  angle=0.0,
                  zorder=100,
                  projection='rect',
                  useRadians=False,
                  org=0,
                  **kwargs):

    verts = regPatch.get_verts()

    if projection != 'rect':
        RA = verts[:, 0]
        RA = np.remainder(RA + 360 - org, 360)  # shift RA values
        ind = RA > 180.
        RA[ind] -= 360  # scale conversion to [-180, 180]
        RA = -RA  # reverse the scale: East to the left

        verts[:, 0] = RA

    if useRadians:
        verts = verts * np.pi / 180.

    ax.plot(verts[:, 0], verts[:, 1], ls='-', c=kwargs['color'], alpha=1, zorder=zorder)

    return regPatch


def plotHSC(ax, **kwargs):

    color = kwargs.get('color', defaultColours['HSC'])
    projection = kwargs.get('projection', 'rect')

    if projection == 'rect':
        regsToPlot = HSC
    else:
        regsToPlot = [HSC_SGC_Mollweide, HSC_3, HSC_4]

    for region in regsToPlot:
        plotPatch(ax, region, color=color, **kwargs)

    addText(ax, 145, -8, 'HSC', color=color, **kwargs)


def plotUKIDSS(ax, **kwargs):

    color = kwargs.get('color', defaultColours['UKIDSS'])
    projection = kwargs.get('projection', 'rect')

    if projection == 'rect':
        regsToPlot = UKIDSS
    else:
        regsToPlot = [UKIDSS_SGC_Mollweide, UKIDSS[2], UKIDSS[3], UKIDSS[4], UKIDSS[5]]

    for region in regsToPlot:
        plotPatch(ax, region, color=color, **kwargs)

    addText(ax, 108, 20, 'UKIDSS', color=color, **kwargs)

    return


def plotApertifMedDeep(ax, **kwargs):

    color = kwargs.get('color', defaultColours['Apertif'])

    # plotRectangle(ax, Perseus-Pisces, **kwargs)
    plotPatch(ax, CVn, color=color, **kwargs)
    plotPatch(ax, HETDEX, color=color, **kwargs)

    addText(ax, 210, 59, 'Apertif', color=color, **kwargs)

    return


def plotALFALFA(ax, **kwargs):

    color = kwargs.get('color', defaultColours['ALFALFA'])
    projection = kwargs.get('projection', 'rect')

    if projection == 'rect':
        regsToPlot = ALFALFA
    else:
        regsToPlot = [ALFALFA_1, ALFALFA_Mollweide]

    for region in regsToPlot:
        plotPatch(ax, region, color=color, **kwargs)

    addText(ax, 110, 31, 'ALFALFA', color=color, **kwargs)

    return


def plotALFALFA_Detailed(ax, **kwargs):

    color = kwargs.get('color', defaultColours['ALFALFA'])
    projection = kwargs.get('projection', 'rect')

    if projection == 'rect':
        regsToPlot = ALFALFA_Detailed
    else:
        regsToPlot = [ALFALFA_1, ALFALFA_Mollweide]

    for region in regsToPlot:
        plotPatch(ax, region, color=color, **kwargs)

    addText(ax, 100, 20, 'ALFALFA', color=color, **kwargs)

    return


def plotATLAS(ax, **kwargs):

    color = kwargs.get('color', defaultColours['ATLAS'])

    plotRectangle(ax, ATLAS, color=color, **kwargs)
    addText(ax, 235, 18, 'H-ATLAS', he='right', color=color, **kwargs)

    return


def plotGAMA(ax, **kwargs):

    color = kwargs.get('color', defaultColours['GAMA'])

    for region in GAMA:
        plotPatch(ax, region, color=color, **kwargs)

    addText(ax, 200, -9, 'GAMA (SAMI)', ha='center', color=color, **kwargs)


# HSC regions
HSC_1 = getRectangle((22 * 15, 360, -1, 7))
HSC_2 = getPolygon(
    np.array([[0, -1], [1.83333 * 15, -1], [1.83333 * 15, -7], [2.66666 * 15, -7],
              [2.66666 * 15, 7], [0, 7], [0, -1]]))

HSC_SGC_Mollweide = getPolygon(
    np.array([[22 * 15, -1], [1.83333 * 15, -1], [1.83333 * 15, -7], [2.66666 * 15, -7],
              [2.66666 * 15, 7], [22 * 15, 7], [22 * 15, -1]]))

HSC_3 = getRectangle((8.5 * 15, 15 * 15, -2, 5))
HSC_4 = getRectangle((13.3 * 15, 16.6666 * 15, 42.5, 44))
HSC = [HSC_1, HSC_2, HSC_3, HSC_4]

HSC_S_Wide_vertices = np.array([[127.5, -2.5], [127.5, 5.5], [170., 5.5], [170., 4.0], [210., 4.0],
                                [210., 5.5], [225., 5.5], [225., -2.5], [127.5, -2.5]])

HSC_S_Wide = getPolygon(HSC_S_Wide_vertices)

# A special region which is the HSC-N regions a bit wider
HSC_N_Wide = getRectangle((15. * 15, 16.6666 * 15, 41.5, 45))

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
ALFALFA_1 = getPolygon(np.array([[110, 0.], [250, 0], [250, 18], [232, 18], [232, 24],
                                 [250, 24], [250, 32], [232, 32], [232, 37], [140, 37],
                                 [140, 32], [110, 32], [110, 24], [140, 24], [140, 18],
                                 [110, 18], [110, 0.]]))
ALFALFA_2 = getRectangle((0, 3. * 15, 0, 36))
ALFALFA_3 = getRectangle((22 * 15, 24 * 15, 0, 36))
ALFALFA = [ALFALFA_1, ALFALFA_2, ALFALFA_3]

ALFALFA_Mollweide = getRectangle((22 * 15, 3 * 15, 0, 36))

# A more detailed footprint of the SGC
ALFALFA_SGC_Detailed_1 = getPolygon(np.array([[324, 36.], [360, 36], [360, 0], [328, 0],
                                              [328, 10], [324, 10], [324, 36]]))
ALFALFA_SGC_Detailed_2 = getPolygon(np.array([[0, 36.], [48, 36], [48., 30], [46, 30],
                                              [46, 22], [48, 22], [48, 14], [38, 14],
                                              [38, 10], [48, 10], [48, 0], [0, 0], [0, 36]]))
ALFALFA_Detailed = [ALFALFA_1, ALFALFA_SGC_Detailed_1, ALFALFA_SGC_Detailed_2]

# HETDEX field
HETDEX_vertices = np.array([[164.0, 45.5], [180.0, 48.9], [196.0, 49.7], [210.0, 48.8],
                            [226.0, 45.0], [229.0, 51.0], [210.0, 55.4], [196.0, 56.2],
                            [180.0, 55.3], [160.5, 51.0], [164.0, 45.5]])

CVn_vertices = np.array([[184.98, 31.55], [191.76, 31.55], [192.51, 46.30], [184.19, 46.30],
                         [184.98, 31.55]])

HETDEX = getPolygon(HETDEX_vertices)
PerseusPisces = getRectangle((23.7, 33.3, 29.9, 37.9))
CVn = getPolygon(CVn_vertices)

# Apertif medium-depth fields
ApertifMedDeep = [HETDEX, PerseusPisces, CVn]

# Latest UKIDSS footprint

UKIDSS_Region1 = np.array([[061.95838, +00.87955], [054.13354, +06.28934], [48.88450, +06.95682],
                           [046.32311, +08.43322], [039.96898, +08.81795], [035.20637, +13.01168],
                           [035.12406, +14.31154], [030.93135, +17.09270], [0.0000000, +17.00000],
                           [0.0000000, -02.00000], [027.91088, -02.19216], [047.03013, -02.55489],
                           [059.19465, -02.59425], [062.08225, -00.05170], [061.95838, +00.87955]])

UKIDSS_Region2 = np.array([[360.00000, +17.00000], [351.87501, +16.96904], [349.93036, +17.04856],
                           [345.16645, +15.03193], [339.41813, +11.29414], [332.81999, +05.71610],
                           [326.77426, +02.16565], [308.18905, +02.31217], [308.04936, -02.10125],
                           [333.13159, -02.30051], [359.89682, -02.32808], [360.00000, +17.00000]])

UKIDSS_Region3 = np.array(
    [[239.36011, -02.72270], [239.49965, +02.40468], [240.41374, +09.81649],
     [237.14271, +11.35246], [220.12234, +15.07941], [209.19142, +15.81317],
     [206.64129, +17.20587], [180.21959, +16.95451], [178.98269, +17.89222],
     [167.77638, +17.63917], [166.22939, +16.14737], [155.82172, +15.56455],
     [148.06986, +14.87760], [123.50277, +09.32339], [124.93499, -02.63832],
     [133.96699, -04.23402], [169.63894, -03.19124], [171.64508, -05.06182],
     [186.01105, -04.97972], [205.37750, -04.99526], [207.43854, -03.45983],
     [225.76219, -03.65258], [232.18386, -03.77643], [239.36011, -02.72270]])

UKIDSS_Region4 = np.array([[211.18558, +21.91363], [211.38403, +37.05334], [188.58584, +36.73161],
                           [189.76090, +21.58410], [198.23781, +20.75589], [211.18558, +21.91363]])

UKIDSS_Region5 = np.array([[242.50593, +20.93451], [230.69122, +27.35972], [220.47165, +32.45357],
                           [222.59255, +35.00932], [226.67155, +34.71249], [232.46691, +32.07994],
                           [234.40908, +31.63117], [238.06585, +28.74299], [240.20232, +31.69502],
                           [244.57836, +32.74295], [252.16101, +33.78028], [255.60497, +30.72946],
                           [250.29212, +25.29416], [245.53296, +21.25760], [242.50593, +20.93451]])

UKIDSS_Region6 = np.array([[120.38225, +17.00487], [114.65075, +20.61480], [110.51386, +25.41934],
                           [112.14745, +27.20383], [111.99926, +28.35569], [114.62997, +30.68592],
                           [132.00404, +30.98818], [135.50603, +33.43509], [136.85960, +34.49213],
                           [143.02820, +36.78138], [147.04215, +36.75095], [149.45191, +34.14314],
                           [147.78351, +31.55026], [145.65717, +29.22833], [138.61797, +26.43600],
                           [130.26806, +22.78999], [125.86197, +19.94148], [121.75807, +17.33957],
                           [120.38225, +17.00487]])

UKIDSS = [getPolygon(UKIDSS_Region1),
          getPolygon(UKIDSS_Region2),
          getPolygon(UKIDSS_Region3),
          getPolygon(UKIDSS_Region4),
          getPolygon(UKIDSS_Region5),
          getPolygon(UKIDSS_Region6)]

# Joined regs 1 and 2 for Mollweide projection
UKIDSS_SGC_Mollweide_vertices = np.array(
    [[0.0000000, -02.32808], [027.91088, -02.19216], [047.03013, -02.55489],
     [059.19465, -02.59425], [062.08225, -00.05170], [061.95838, +00.87955],
     [054.13354, +06.28934], [48.88450, +06.95682], [046.32311, +08.43322],
     [039.96898, +08.81795], [035.20637, +13.01168], [035.12406, +14.31154],
     [030.93135, +17.09270], [360.00000, +17.00000], [351.87501, +16.96904],
     [349.93036, +17.04856], [345.16645, +15.03193], [339.41813, +11.29414],
     [332.81999, +05.71610], [326.77426, +02.16565], [308.18905, +02.31217],
     [308.04936, -02.10125], [333.13159, -02.30051], [359.89682, -02.32808],
     [360.00000, -02.32808]])

UKIDSS_SGC_Mollweide = getPolygon(UKIDSS_SGC_Mollweide_vertices)

# Some renaming, for convenience
UKIDSS_240 = UKIDSS[4]
UKIDSS_120 = UKIDSS[5]
UKIDSS_ATLAS = UKIDSS[3]
ALFALFA_NGC = ALFALFA[0]


def getPlatesInFootprint(plates, coords=False):
    """Returns the list of plates that are within MaNGA footprint."""

    from Totoro.dbclasses.plate import Plate

    plates = np.atleast_1d(plates)

    if coords:
        tmpPlates = [Plate.createMockPlate(ra=plate[0], dec=plate[1]) for plate in plates]
        plates = tmpPlates

    footprintPlates = []
    for plate in plates:

        # Excludes the two most eastern GAMA fields
        if (GAMA_2.get_path().contains_points([plate.coords])[0] or
                GAMA_3.get_path().contains_points([plate.coords])[0]):
            continue

        if ((UKIDSS_120.get_path().contains_points([plate.coords])[0] and
             ALFALFA_NGC.get_path().contains_points([plate.coords])[0]) or
            (UKIDSS_240.get_path().contains_points([plate.coords])[0] and
             ALFALFA_NGC.get_path().contains_points([plate.coords])[0]) or
                ATLAS.get_path().contains_points([plate.coords])[0] or
                SGC[0].get_path().contains_points([plate.coords])[0] or
                SGC[1].get_path().contains_points([plate.coords])[0] or
                HETDEX.get_path().contains_points([plate.coords])[0] or
                PerseusPisces.get_path().contains_points([plate.coords])[0] or
                CVn.get_path().contains_points([plate.coords])[0] or
                HSC_S_Wide.get_path().contains_points([plate.coords])[0] or
                HSC_N_Wide.get_path().contains_points([plate.coords])[0]):

            footprintPlates.append(plate)

    if len(footprintPlates) == 1:
        return footprintPlates[0]
    else:
        return footprintPlates


def plotFootprint(ax, regions='all', projection='rect', org=0):
    """Overplots the footprint on a Matplotlib axes instance.

    Parameters
    ----------
    ax : matplotlib axes
        The axes on which to overplot the regions.

    regions : string, None, list of strings
        The list of footprint regions to overplot. If `'all'` or `None`, all
        the regions are plotted. A list of regions can also be provided.

    proj : string
        The projection to use, either `'rect'` or `'mollweide'`

    org : float
        The origin of the projection.

    Returns
    -------
    ax : matplotlib axes
        The same input axes, `ax`, after the regions have been overplot.

    """

    if not regions or regions == 'all':
        plots = [plotGAMA, plotApertifMedDeep, plotATLAS, plotHSC, plotUKIDSS, plotALFALFA]
    else:
        if isinstance(regions, str):
            regions = [regions]
        plots = []
        for reg in regions:
            assert 'plot{0}'.format(reg) in globals(), '{0} is not defined'.format(reg)
            plots.append(eval('plot{0}'.format(reg)))

    useRadians = True if projection != 'rect' else False

    for plot in plots:
        plot(ax, projection=projection, useRadians=useRadians, org=org)

    return ax


# def plotEBHIS(ax, **kwargs):
#
#     ax.axhline(-6, color=kwargs['color'], lw=kwargs['lw'], ls='dashed',
#                zorder=10)
#     ax.annotate('', xy=(260, -5.5), xytext=(260, 6), textcoords='data',
#                 arrowprops=dict(facecolor=kwargs['color'],
#                                 edgecolor=kwargs['color'],
#                                 linewidth=kwargs['lw'],
#                                 arrowstyle='<-'), zorder=20)
#     ax.text(262, -2, 'EBHIS', color=kwargs['color'],
#             fontsize=12, zorder=20, weight='heavy',
#             fontdict={'family': 'sans-serif'})
#
#     return
#
#
# def plotApertifShallow(ax, **kwargs):
#
#     ax.axhline(27, color=kwargs['color'], lw=kwargs['lw'], ls='dashed',
#                zorder=10)
#     ax.annotate('', xy=(300, 27.5), xytext=(300, 37.5), textcoords='data',
#                 arrowprops=dict(facecolor=kwargs['color'],
#                                 edgecolor=kwargs['color'],
#                                 linewidth=kwargs['lw'],
#                                 arrowstyle='<-'), zorder=20)
#     ax.text(260, 39.5, 'Apertif (shallow)', color=kwargs['color'],
#             fontsize=12, ha='left', zorder=20, weight='heavy',
#             fontdict={'family': 'sans-serif'})
#
#     return
#
#
# def plotASKAP(ax, **kwargs):
#
#     ax.axhline(30, color=kwargs['color'], lw=kwargs['lw'], ls='dashed',
#                zorder=10)
#     ax.annotate('', xy=(80, 29.5), xytext=(80, 14.5), textcoords='data',
#                 arrowprops=dict(facecolor=kwargs['color'],
#                                 edgecolor=kwargs['color'],
#                                 linewidth=kwargs['lw'],
#                                 arrowstyle='<-'), zorder=20)
#     ax.text(99, 11, 'ASKAP', color=kwargs['color'], fontsize=12, ha='right',
#             zorder=20, weight='heavy', fontdict={'family': 'sans-serif'})
#
#     return
