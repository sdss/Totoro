#!/usr/bin/env python
# encoding: utf-8
"""
createSchedule.py

Created by José Sánchez-Gallego on 8 Sep 2014.
Licensed under a 3-clause BSD license.

Revision history:
    8 Sep 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import sys
from astropy import table
from astropy import time


def getChange(change):

    if change == '-999.':
        raise ValueError('change=-999.')

    change = map(float, change.split(':'))
    change = change[0] + change[1] / 60.

    return time.TimeDelta(change * 3600, format='sec')


def createSchedule(finalSchedule):

    finalScheduleTable = table.Table.read(
        finalSchedule, format='ascii.commented_header')

    formattedTable = table.Table(
        None, names=['JD', 'ID1', 'ID2', 'Fraction', 'Value1', 'Value2',
                     'ApogeeStart', 'ApogeeEnd', 'eBossStart', 'eBossEnd',
                     'mangaStart', 'mangaEnd', 'ID3'],
        dtype=[float, int, int, float, float, float, float, float,
               float, float, float, float, int])

    for night in finalScheduleTable:

        mjd = night['MJD']

        eveTwUT = map(float, night['eveTwUT'].split(':'))
        eveTwUT = eveTwUT[0] + eveTwUT[1] / 60.

        morTwUT = map(float, night['morTwUT'].split(':'))
        morTwUT = morTwUT[0] + morTwUT[1] / 60.

        changeUT = night['ChangeUT']

        tt = time.Time(mjd, format='mjd')
        eveTw = tt + time.TimeDelta(eveTwUT * 3600, format='sec')
        morTw = tt + time.TimeDelta(morTwUT * 3600, format='sec')

        if night['orApogee2'] == 1:
            mangaStart = 0
            mangaEnd = 0
            eBossStart = 0
            eBossEnd = 0
            if night['orEboss'] == 0 and night['orManga'] == 0:
                apogeeStart = eveTw.jd
                apogeeEnd = morTw.jd
            elif night['orEboss'] == 2:
                change = tt + getChange(changeUT)
                apogeeStart = eveTw.jd
                apogeeEnd = change.jd
                eBossStart = change.jd
                eBossEnd = morTw.jd
            elif night['orManga'] == 2:
                change = tt + getChange(changeUT)
                apogeeStart = eveTw.jd
                apogeeEnd = change.jd
                mangaStart = change.jd
                mangaEnd = morTw.jd
        elif night['orEboss'] == 1:
            change = tt + getChange(changeUT)
            eBossStart = eveTw.jd
            eBossEnd = change.jd
            if night['orManga'] == 2:
                mangaStart = change.jd
                mangaEnd = morTw.jd
                apogeeStart = 0
                apogeeEnd = 0
            elif night['orApogee2'] == 2:
                apogeeStart = change.jd
                apogeeEnd = morTw.jd
                mangaStart = 0
                mangaEnd = 0
        elif night['orManga'] == 1:
            change = tt + getChange(changeUT)
            mangaStart = eveTw.jd
            mangaEnd = change.jd
            if night['orEboss'] == 2:
                eBossStart = change.jd
                eBossEnd = morTw.jd
                apogeeStart = 0
                apogeeEnd = 0
            elif night['orApogee2'] == 2:
                apogeeStart = change.jd
                apogeeEnd = morTw.jd
                eBossStart = 0.
                eBossEnd = 0.

        newRow = (int(tt.jd), 0, 0, 0.0, 0.0, 0.0, apogeeStart, apogeeEnd,
                  eBossStart, eBossEnd, mangaStart, mangaEnd, 2)

        formattedTable.add_row(newRow)

    formattedTable.write('tmpSchedule.dat',
                         format='ascii.fixed_width_no_header', delimiter=' ')


if __name__ == '__main__':
    createSchedule(sys.argv[1])
