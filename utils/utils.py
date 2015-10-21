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
from sdss.internal.manga.Totoro import log
from sdss.internal.manga.Totoro import config
from sdss.internal.manga.Totoro import exceptions
from sdss.internal.manga.Totoro.apoDB import TotoroDBConnection
from sdss.manga.mlhalimit import mlhalimit as mlhalimitHours
from collections import OrderedDict
from itertools import combinations
import numpy as np
import warnings
import os


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


def isPlateComplete(plate, format='plate_id', forceCheckCompletion=False,
                    **kwargs):
    """Returns if a plate is complete using the MaNGA logic. If
    forceCheckCompletion is False and the plugging is marked as complete,
    no plateCompletion check is performed (this saves some time)."""

    from sdss.internal.manga.Totoro.dbclasses.plate import Plate

    if not isinstance(plate, Plate):
        if format.lower() not in ['pk', 'plate_id']:
            raise exceptions.TotoroError('format must be plate_id or pk.')
        plate = Plate(plate, format=format.lower(), **kwargs)

    plugStatus = [plugging.status.label for plugging in plate.pluggings]

    if 'Good' in plugStatus or 'Overridden Good' in plugStatus:
        plugComplete = True
    elif u'Overridden Incomplete' in plugStatus:
        plugComplete = False
    else:
        plugComplete = None

    if plugComplete is not None and forceCheckCompletion is False:
        return plugComplete

    if plate.getPlateCompletion(includeIncompleteSets=False) >= 1.:
        plateComplete = True
    else:
        plateComplete = False

    if (plateComplete is True and
            np.isnan(plate.getCumulatedSN2(includeIncomplete=False)).any()):
        log.debug('plate_id={0}: not all cameras have been correctly reduced. '
                  'Setting plateComplete=False.'.format(plate.plate_id))
        plateComplete = False

    if plugComplete is not None:
        if plugComplete is not plateComplete:
            warnings.warn('plate={0}: plugging status is {1} but calculated '
                          'status is {2}.'.format(
                              plate.plate_id,
                              'complete' if plugComplete else 'incomplete',
                              'complete' if plateComplete else 'incomplete'),
                          exceptions.TotoroUserWarning)
            return plugComplete
        else:
            return plugComplete
    else:
        return plateComplete


def getAPOcomplete(plates, format='plate_id',
                   SN2_blue=None, SN2_red=None, limitSN=False,
                   func=np.max, createFile=False, **kwargs):
    """Returns a dictionary with the APOcomplete output."""

    from sdss.internal.manga.Totoro.dbclasses import Plate

    SN2_blue = config['SN2thresholds']['plateBlue'] \
        if SN2_blue is None else SN2_blue
    SN2_red = config['SN2thresholds']['plateRed'] \
        if SN2_red is None else SN2_red

    format = format.lower()
    if format.lower() not in ['pk', 'plate_id']:
        raise exceptions.TotoroError('format must be plate_id or pk.')

    plates = np.atleast_1d(plates)

    APOcomplete = OrderedDict()

    for plate in plates:

        setsToAPOcomplete = None

        if not isinstance(plate, Plate):
            plate = Plate(plate, format=format.lower(), **kwargs)

        if isPlateComplete(plate) is False:
            warnings.warn('plate_id={0} is not complete. APOcomplete output '
                          'must not be used.'.format(plate.plate_id),
                          exceptions.TotoroUserWarning)

        APOcomplete[plate.plate_id] = []

        validSets = plate.getValidSets()

        if limitSN:
            for nSets in range(1, len(validSets)+1):

                combSets = list(combinations(validSets, nSets))
                SN2 = np.array([_cumulatedSN2(sets) for sets in combSets])

                overSN2 = [(combSets[ii], SN2[ii])
                           for ii in range(len(combSets))
                           if SN2[ii][0] >= SN2_blue and SN2[ii][1] >= SN2_red]

                if len(overSN2) > 0:

                    relativeSN2 = np.array(
                        [(overSN2[ii][1][0] / SN2_blue) *
                         (overSN2[ii][1][1] / SN2_red)
                         for ii in range(len(overSN2))])

                    setsToAPOcomplete = overSN2[
                        np.where(relativeSN2 == func(relativeSN2))[0][0]][0]

                    break

        if setsToAPOcomplete is None:
            setsToAPOcomplete = validSets

        for ss in setsToAPOcomplete:
            for exp in ss.totoroExposures:

                mjd = exp.getMJD()
                pk = ss.pk
                dPos = exp.ditherPosition.upper()
                nExp = exp.exposure_no

                APOcomplete[plate.plate_id].append(
                    [plate.plate_id, mjd, pk, dPos, nExp])

        if len(setsToAPOcomplete) > 0:
            apoCompleteSN2 = _cumulatedSN2(setsToAPOcomplete)
        else:
            apoCompleteSN2 = np.array([0., 0.])

        if apoCompleteSN2[0] >= SN2_blue and apoCompleteSN2[1] >= SN2_red:
            log.info('APOcomplete for plate_id={0} generated with '
                     'SN2_blue={1:.1f}, SN2_red={2:.1f}.'.format(
                         plate.plate_id, apoCompleteSN2[0], apoCompleteSN2[1]))
        else:
            warnings.warn('APOcomplete for plate_id={0} generated with '
                          'SN2_blue={1:.1f}, SN2_red={2:.1f}, which is lower '
                          'than the thresholds.'.format(
                              plate.plate_id, apoCompleteSN2[0],
                              apoCompleteSN2[1]), exceptions.TotoroUserWarning)

    if createFile:
        createAPOcompleteFile(APOcomplete, **kwargs)

    return APOcomplete


def createAPOcompleteFile(APOcomplete, path=None):
    """Writes the APOcomplete file in Yanny format."""

    path = './' if path is None else path

    for key in APOcomplete:

        apocompPath = os.path.join(path, 'apocomp-{0:04d}.par'.format(key))

        data = APOcomplete[key]

        plateid = [dd[0] for dd in data]
        mjd = [dd[1] for dd in data]
        setno = [dd[2] for dd in data]
        dpos = [dd[3] for dd in data]
        expno = [dd[4] for dd in data]

        strstruct = 'typedef struct {\n long plateid;\n long mjd;\n ' + \
            'int set;\n char mgdpos[2];\n long exposure;\n} APOCOMP;\n\n'

        ff = open(apocompPath, 'w')
        ff.write(strstruct)

        for ii in range(len(data)):
            expstr = '{0} {1} {2} {3} {4} {5}\n'.format(
                'APOCOMP', plateid[ii], mjd[ii], setno[ii],
                dpos[ii], expno[ii])
            ff.write(expstr)

        ff.close()

    return


def _cumulatedSN2(sets):
    """Returns the cumulated SN2 for a list of sets as [SN2_blue, SN2_red]."""

    SN2array = np.array([ss.getSN2Array() for ss in sets])
    if len(SN2array) == 0:
        SN2array = np.array([[0.0, 0.0, 0.0, 0.0]])

    SN2sum = np.sum(SN2array, axis=0)

    return np.array([np.mean(SN2sum[0:2]), np.mean(SN2sum[2:])])


def JDdiff(JD0, JD1):
    """Returns the number of seconds between two Julian dates."""

    return (JD1 - JD0) * 86400


def isMaNGA_Led(plate):
    """Returns True if the plate is a MaNGA-led plate."""

    totoroDB = TotoroDBConnection()
    plateDB = totoroDB.plateDB
    session = totoroDB.session

    from sdss.internal.manga.Totoro.dbclasses import Plate

    if isinstance(plate, (Plate, plateDB.Plate)):
        pass
    else:
        try:
            with session.begin(subtransactions=True):
                plate = session.query(plateDB.Plate).filter(
                    plateDB.Plate.plate_id == plate).one()
        except:
            return False

    for survey in plate.surveys:
        if (survey.label == 'MaNGA' and plate.currentSurveyMode is not None and
                plate.currentSurveyMode.label == 'MaNGA dither'):
            return True

    return False
