#!/usr/bin/env python
# encoding: utf-8
"""
conversion.py

Created by José Sánchez-Gallego on 13 Dec 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

This file includes classes and function to convert
between different types of data.

"""

__ALL__ = ['jd2lmst', 'lmst2time']

from astropy import time
from astropy import coordinates as coo
from astropy import units as uu
from ..core.defaults import LONGITUDE


def jd2lmst(jd, longitude=LONGITUDE()):

    if isinstance(jd, time.Time):
        jd = jd.jd

    jd0 = int(jd) + 0.5
    dd0 = jd0 - 2451545.
    dd = jd - 2451545.
    tt = dd / 36525
    hh = (jd - jd0) * 24.

    gmst = 6.697374558 + 0.06570982441908 * dd0 + 1.00273790935 * hh + \
        0.000026 * tt**2
    gmstL = coo.Longitude(gmst * uu.hour)

    if not isinstance(longitude, coo.Longitude):
        longitude = coo.Longitude(longitude * uu.degree)

    return coo.Longitude(gmstL + longitude)


def lmst2time(lmst, mjd, longitude=LONGITUDE()):
    """Inverse function for utc2lmst().

    Parameters
    ----------
    lmst : float
        The LMST in hours.
    mjd : float
        The modified julian day of the observation.
    longitude : float or `astropy.coordinates.angles.Longitude`, optional
        The longitude of the site at which the LMST
        is calculated. If not define, APO longitude will
        be used. East longitudes with longitude in the range
        [0, 360) degrees should be used.

    Returns
    -------
    utc : `astropy.time.Time` object
        A time object with the time corresponding to the pair mjd and lmst.

    """

    lmst = coo.Longitude(lmst, unit=uu.hour)

    if not isinstance(longitude, coo.Longitude):
        longitude = coo.Longitude(longitude, unit=uu.degree)

    mjd = time.Time(mjd, format='mjd', scale='utc')

    # Calculates the sidereal time at the previous midnight
    midNight = time.Time(mjd.jd1, format='jd', scale='utc')
    midNightSidTime = jd2lmst(midNight, longitude=longitude)

    # Calculates UTC
    utc = lmst - midNightSidTime

    # Corrects for the difference between sidereal and UTC seconds
    utc = utc * 365.25 / 366.25

    # Converts UTC to a TimeDelta
    utc = time.TimeDelta(utc.hour * 3600, format='sec')

    # Creates a Time object from the MJD and UTC
    utcDate = midNight + utc

    return utcDate
