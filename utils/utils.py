#!/usr/bin/env python
# encoding: utf-8
"""
utils.py

Created by José Sánchez-Gallego on 15 May 2014.
Licensed under a 3-clause BSD license.

Revision history:
    15 May 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import numpy as np
from .. import config
from astropysics import obstools
from .. import log
from ..exceptions import TotoroError, TotoroUserWarning
import warnings
from sdss.manga import mlhalimit as mlhalimitHours
from collections import OrderedDict
from astropy import time, units


def mlhalimit(dec):
    """Returns HA limits in degrees."""

    return mlhalimitHours(dec) * 15.


def computeAirmass(dec, ha, lat=config['observatory']['latitude'],
                   correct=[75., 10.]):
    """Calculates the airmass for a given declination and HA (in degrees).

    By default, assumes that the latitude of the observation is the one set
    in the configuration file. If correct is defined, abs(HA) anggles greater
    than correct[0] are given a flat value correct[1].
    """

    dec = np.atleast_1d(dec)
    ha = np.atleast_1d(ha) % 360.

    if ha > 180:
        ha -= 360

    airmass = (np.sin(lat * np.pi / 180.) * np.sin(dec * np.pi / 180.) +
               np.cos(lat * np.pi / 180.) * np.cos(dec * np.pi / 180.) *
               np.cos(ha * np.pi / 180.)) ** (-1)

    if correct is not None:
        airmass[np.abs(ha) > correct[0]] = correct[1]

    if len(airmass) == 1:
        return airmass[0]
    else:
        return airmass


def isPlateComplete(plate, format='plate_id', **kwargs):
    """Returns if a plate is complete using the MaNGA logic."""

    from ..dbclasses import Plate

    if not isinstance(plate, Plate):
        if format.lower() not in ['pk', 'plate_id']:
            raise TotoroError('format must be plate_id or pk.')
        plate = Plate(plate, format=format.lower(), **kwargs)

    if plate.getPlateCompletion(includeIncompleteSets=False) > 1.:
        plateComplete = True
    else:
        plateComplete = False

    plugStatus = np.array(
        [plugging.status.label for plugging in plate.pluggings])
    if 'Good' in plugStatus or 'Overriden Good' in plugStatus:
        plugComplete = True
    elif 'Overriden Incomplete' in plugStatus:
        plugComplete = False
    else:
        plugComplete = None

    if plugComplete is not None:
        if plugComplete is not plateComplete:
            warnings.warn('plugging status is {0} but calculated '
                          'status is {1}.'.format(
                              'complete' if plugComplete else 'incomplete',
                              'complete' if plateComplete else 'incomplete'),
                          TotoroUserWarning)
            return plugComplete
    else:
        return plateComplete


def getAPOcomplete(plates, format='plate_id', **kwargs):
    """Returns a dictionary with the APOcomplete output."""

    from ..dbclasses import Plate

    format = format.lower()
    if format.lower() not in ['pk', 'plate_id']:
        raise TotoroError('format must be plate_id or pk.')

    plates = np.atleast_1d(plates)

    APOcomplete = OrderedDict()

    for plate in plates:

        if not isinstance(plate, Plate):
            plate = Plate(plate, format=format.lower(), **kwargs)

        if isPlateComplete(plate) is False:
            warnings.warn('plate_id={0} is not complete. APOcomplete output '
                          'must not be used.'.format(plate.plate_id),
                          TotoroUserWarning)

        APOcomplete[plate.plate_id] = []

        for set in plate.sets:
            for exp in set.totoroExposures:

                mjd = exp.getMJD()
                ss = set.pk
                dPos = exp.ditherPosition.upper()
                nExp = exp.exposure_no

                APOcomplete[plate.plate_id].append(
                    [plate.plate_id, mjd, ss, dPos, nExp])

    return APOcomplete


def createSite(longitude=None, latitude=None, altitude=None,
               name=None, verbose=False, **kwargs):
    """Returns an astropysics.obstools.site instance. By default, uses the
    coordinates for APO."""

    if None in [longitude, latitude, altitude, name]:
        assert 'observatory' in config.keys()

    longitude = config['observatory']['longitude'] \
        if longitude is None else longitude
    latitude = config['observatory']['latitude'] \
        if latitude is None else latitude
    altitude = config['observatory']['altitude'] \
        if altitude is None else altitude

    if name is None:
        if 'name' not in config['observatory']:
            name = ''
        else:
            name = config['observatory']['name']

    site = obstools.Site(latitude, longitude, name=name, alt=altitude)
    # site.localSiderialTime = lambda jd: _calculateLST(site, jd)

    if verbose:
        log.info('Created site with name \'{0}\''.format(name))

    return site


# def _calculateLST(site, jd):
#     """Not-to-be-used-directly function to replace the imprecise
#     locaSiderialTime in astropysics."""

#     tmpTime = time.Time(jd, format='jd', scale='tai')
#     tmpTime.delta_ut1_utc = 0.

#     lst = tmpTime.sidereal_time('apparent',
#                                 longitude=float(site.longitude.degrees))

#     return lst.hour


def JDdiff(JD0, JD1):
    """Returns the number of seconds between two Julian dates."""

    return (JD1 - JD0) * 86400
