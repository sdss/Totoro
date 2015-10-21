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

from astropy.coordinates.angles import Longitude
from ..exceptions import TotoroError
from ..core import ConfigObject
from astropy.units import degree, hour, hourangle
from astropy.coordinates import Angle
from astropy import time

APO_LONGITUDE = 254.179722
LONGITUDE = ConfigObject('longitude', APO_LONGITUDE,
                         'The longitude of the observatory')


def utc2lmst(utc, longitude=LONGITUDE()):
    """Returns the LMST for a UTC time.

    This function returns the LMST for a certain UTC time
    assuming that delta_ut1_utc=0.0.

    Parameters
    ----------
    utc : `astropy.time.Time` object
        The UTC time, in astropy format.
    longitude : float or `astropy.coordinates.angles.Longitude` object
        The longitude of the site at which the LMST
        is calculated. If not define, APO longitude will
        be used. East longitudes with longitude in the range
        [0, 360) degrees should be used.

    Returns
    -------
    lmst : `astropy.coordinates.angle.Longitude` object
        The LMST for the input UTC at the specified longitude.

    """

    if isinstance(longitude, Longitude):
        pass
    else:
        try:
            longitude = Longitude(longitude, unit=degree)
        except:
            raise TotoroError('longitude cannot be understood.')

    utc.delta_ut1_utc = 0.0
    utc.lon = longitude

    return utc.sidereal_time('mean')


def lmst2utc(lmst, mjd, longitude=LONGITUDE()):
    """Inverse function for utc2lmst().

    Parameters
    ----------
    lmst : `astropy.coordinates.angle.Longitude` object or float
        The LMST in astropy Longitude format or as a float in hours.
    longitude : float or `astropy.coordinates.angles.Longitude` object
        The longitude of the site at which the LMST
        is calculated. If not define, APO longitude will
        be used. East longitudes with longitude in the range
        [0, 360) degrees should be used.

    Returns
    -------
    utc : `astropy.time.Time` object
        A time object with the time corresponding to the pair mjd and lmst.

    """

    if not isinstance(lmst, Longitude):
        lmst = Angle(lmst, unit=hourangle)

    if isinstance(longitude, Longitude):
        pass
    else:
        try:
            longitude = Longitude(longitude, unit=degree)
        except:
            raise TotoroError('longitude cannot be understood.')

    # Rewraps longitude to the range -180 to 180
    longitude.wrap_angle = 180 * degree

    if isinstance(mjd, time.Time):
        pass
    else:
        mjd = time.Time(mjd, format='mjd', scale='utc')

    # Calculates the sidereal time at the previous midnight
    midNight = time.Time(mjd.jd1, format='jd', scale='utc')
    midNightSidTime = utc2lmst(midNight, longitude=longitude)

    # Calculates UTC
    utc = lmst - midNightSidTime

    if utc.hour > 24.:
        utc -= Angle(24, unit=hour)
    elif utc.hour < 0.:
        utc += Angle(24, unit=hour)

    # Corrects for the difference between sidereal and UTC seconds
    utc = utc * 365.25 / 366.25

    # Converts UTC to a TimeDelta
    utc = time.TimeDelta(utc.hour * 3600, format='sec')

    # Creates a Time object from the MJD and UTC
    utcDate = midNight + utc

    return utcDate
