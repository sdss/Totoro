#!/usr/bin/env python
# encoding: utf-8
"""
Site.py

Created by José Sánchez-Gallego on 24 Oct 2014.
Licensed under a 3-clause BSD license.

Revision history:
    24 Oct 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division, print_function

import numpy as np
from astropy import time


# TBD: we need ephem at APO, since our older astropy does not contain
# astropy.coordinates.get_sun()
# TBD: once we update astropy at APO to 1.0, we can remove the ephem dependency.
# for now, we wrap it like this, which will make getSunAltitude break
# if pyephem is not installed.
try:
    import ephem
except Exception:
    ephem = None


class Site(object):
    """A class similar to `astropysics.obstools.Site` to perform LST->calendar
    conversions.

    For now it contains only date->LST and LST->date methods, but any
    conversion tool to be used by Petunia, Totoro and the autoscheduler should
    go here and be used for internal consistency.

    Parameters
    ----------
    longitude : float, optional
        The longitude of the observing site. If None, it uses APO longitude.

    latitude : float, optional
        Ditto for latitude

    altitude : float, optional
        Ditto for the altitude of the site

    name : str, optional
        The name of the observing site

    """

    def __init__(self,
                 longitude=None,
                 latitude=None,
                 altitude=None,
                 name=None,
                 verbose=False,
                 **kwargs):

        self.longitude = 254.179722 if longitude is None else float(longitude)
        self.latitude = 32.766666667 if latitude is None else float(latitude)
        self.altitude = 2788 if altitude is None else float(altitude)

        self.name = 'APO' if name is None else name

        # TBD: we need ephem at APO, since our older astropy does not contain
        # astropy.coordinates.get_sun()
        if ephem:
            # The Sun, as an PyEphem object.
            self.sun = ephem.Sun()
        else:
            self.sun = None

        if verbose:
            print('Created site with name \'{0}\''.format(self.name))

    def localSiderialTime(self, *args, **kwargs):
        """Alias for localSiderealTime observing the wrong spelling in
        astropysics."""

        return self.localSiderealTime(*args, **kwargs)

    def localSiderealTime(self, inputDate=None, format='jd', **kwargs):
        """Returns the LST for a given date.

        Parameters
        ----------
        inputDate : optional
            An `astropy.Time.time` instance or the argument to create one.
            If None, the current time will be used. The UTC scale is used.

        format : string, optional
            If date is not None or a Time instance, the value to be passed
            to the `astropy.time.Time` `format` keyword. The default is `jd`.

        Returns
        -------
        result : float
            A float with the LST time for the given date(s), in hours.

        """

        if inputDate is None:
            inputDate = time.Time.now()

        if isinstance(inputDate, time.Time):
            JD = inputDate.jd
        elif format == 'jd':
            JD = np.array(inputDate)
        else:
            try:
                inputDate = time.Time(inputDate, scale='tai', format=format)
                inputDate.jd
            except Exception:
                raise ValueError('inputDate format not recognised.')

        # inputDate.delta_ut1_utc = 0.
        # lst = inputDate.sidereal_time('apparent', longitude=self.longitude)

        # Days since J2000.
        dd = JD - 2451545.0

        lmst = ((280.46061837 + 360.98564736629 * dd + 0.000388 *
                 (dd / 36525.)**2 + self.longitude) % 360) / 15.

        return lmst

    def localSiderealTimeToDate(self, lst, date=None, format=None):
        """Returns an `astropy.time.Time` instance for a LST
        at a given date.

        Parameters
        ----------
        lst : float or iterable
            A LST value or list of them to be converted to dates.

        date : optional
            An `astropy.Time.time` instance or the argument to create one.
            If None, the current time will be used.

        format : string, optional
            If date is not None or a Time instance, the value to be passed
            to the `astropy.time.Time` `format` keyword.

        Returns
        -------
        result : `astropy.time.Time`
            An `astropy.time.Time` instance of the same size of the lst
            input list, with each element being the date of the corresponding
            LST.

        """

        if date is None:
            date = time.Time.now()

        lst = np.atleast_1d(lst)

        if isinstance(date, time.Time):
            pass
        else:
            try:
                date = time.Time(date, format=format, scale='tai')
            except Exception:
                raise ValueError('date format not recognised.')

        LST0 = self.localSiderealTime(date)
        testPoint = lst[0]
        diffs = (lst - testPoint) % 24

        if np.abs(testPoint - LST0) > 12.:
            delta = (testPoint - LST0) % 24
        else:
            delta = testPoint - LST0

        lstDelta = diffs + delta

        UTDates = date + time.TimeDelta(lstDelta * 3600, format='sec', scale='tai')

        if len(UTDates) == 1:
            return UTDates[0]
        else:
            return UTDates

    def getAltitude(self, timeObject, ra, dec):
        """Returns the altitude in degrees of a `(ra, dec)` target at a
        certain time.

        Parameters
        ----------
        timeObject : `astropy.time.Time` object
            The datetime at which the altitude will be calculated.

        ra, dec : float
            The coordinates of the object.

        Returns
        -------
        result : float
            The altitude of the object at `timeObject` in degrees.

        """

        LST = self.localSiderealTime(timeObject)

        HA = (LST * 15. - ra) % 360.

        sinAlt = (
            np.sin(np.deg2rad(dec)) * np.sin(np.deg2rad(self.latitude)) +
            np.cos(np.deg2rad(dec)) * np.cos(np.deg2rad(self.latitude)) * np.cos(np.deg2rad(HA)))

        return np.rad2deg(np.arcsin(sinAlt))

    def getSunAltitude(self, timeObject=None):
        """Returns the altitude os the Sun at a certain time.

        Parameters
        ----------
        timeObject : `astropy.time.Time` object, optional
            The datetime at which the altitude will be calculated. If not
            defined, the current time will be used.

        Returns
        -------
        result : float
            The altitude of the object at `timeObject` in degrees.

        """

        timeObject = timeObject if timeObject is not None else time.Time.now()

        sunRA = []
        sunDec = []

        isoTimes = np.atleast_1d(timeObject.iso)

        # TBD: we need ephem at APO, since our older astropy does not contain
        # astropy.coordinates.get_sun()
        if not ephem:
            raise RuntimeError('pyephem not installed: cannot compute sun angles.')
        for tt in isoTimes:
            # PyEphem does not accept time arrays.
            self.sun.compute(tt)
            sunRA.append(self.sun.ra * 180. / np.pi)
            sunDec.append(self.sun.dec * 180. / np.pi)

        # Squeezes RA, Dec dimensionality
        sunRA = np.squeeze(sunRA)
        sunDec = np.squeeze(sunDec)

        return self.getAltitude(timeObject, sunRA, sunDec)
