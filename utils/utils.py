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
from sdss.internal.manga.Totoro import config
from sdss.internal.manga.Totoro import exceptions
import warnings
from sdss.manga import mlhalimit as mlhalimitHours
from collections import OrderedDict
from sdss.utilities import Site


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

    from sdss.internal.manga.Totoro.dbclasses import Plate

    if not isinstance(plate, Plate):
        if format.lower() not in ['pk', 'plate_id']:
            raise exceptions.TotoroError('format must be plate_id or pk.')
        plate = Plate(plate, format=format.lower(), **kwargs)

    if plate.getPlateCompletion(includeIncompleteSets=False) > 1.:
        plateComplete = True
    else:
        plateComplete = False

    plugStatus = [plugging.status.label for plugging in plate.pluggings]

    if 'Good' in plugStatus or 'Overridden Good' in plugStatus:
        plugComplete = True
    elif u'Overridden Incomplete' in plugStatus:
        plugComplete = False
    else:
        plugComplete = None

    if plugComplete is not None:
        if plugComplete is not plateComplete:
            warnings.warn('plugging status is {0} but calculated '
                          'status is {1}.'.format(
                              'complete' if plugComplete else 'incomplete',
                              'complete' if plateComplete else 'incomplete'),
                          exceptions.TotoroUserWarning)
            return plugComplete
    else:
        return plateComplete


def getAPOcomplete(plates, format='plate_id', **kwargs):
    """Returns a dictionary with the APOcomplete output."""

    from sdss.internal.manga.Totoro.dbclasses import Plate

    format = format.lower()
    if format.lower() not in ['pk', 'plate_id']:
        raise exceptions.TotoroError('format must be plate_id or pk.')

    plates = np.atleast_1d(plates)

    APOcomplete = OrderedDict()

    for plate in plates:

        if not isinstance(plate, Plate):
            plate = Plate(plate, format=format.lower(), **kwargs)

        if isPlateComplete(plate) is False:
            warnings.warn('plate_id={0} is not complete. APOcomplete output '
                          'must not be used.'.format(plate.plate_id),
                          exceptions.TotoroUserWarning)

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


def JDdiff(JD0, JD1):
    """Returns the number of seconds between two Julian dates."""

    return (JD1 - JD0) * 86400
