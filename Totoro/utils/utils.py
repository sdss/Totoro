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

from __future__ import division, print_function

import os
import subprocess
import warnings
from builtins import range, str
from collections import OrderedDict
from itertools import combinations

import numpy as np
from pydl.pydlutils.yanny import yanny
from scipy.spatial.distance import pdist
from sqlalchemy.exc import InvalidRequestError, ResourceClosedError

from Totoro import config, exceptions, log
from Totoro.db import getConnection


_avoid_cart2_cahche = {}


def mlhalimit(dec):
    """Returns HA limits in DEGREES.

    Calculates the maximum HAs acceptable for a list of declinations.
    Uses the polynomial fit by David Law and a omega limit of 0.5.
    """

    isIterable = hasattr(dec, '__iter__')

    funcFit = np.array([1.78693, 0.0663050, -0.00174096, 2.62002e-05, -1.03959e-07,
                        -1.49150e-09])[::-1]

    dec = np.atleast_1d(dec)
    halimit = np.abs(np.polyval(funcFit, dec)) * 15.

    # halimit[np.where((dec < -10) | (dec > 80))] = 0.0

    return halimit[0] if not isIterable else halimit


def computeAirmass(dec, ha, lat=config['observatory']['latitude'], correct=[75., 10.]):
    """Calculates the airmass for a given declination and HA (in degrees).

    By default, assumes that the latitude of the observation is the one set
    in the configuration file. If correct is defined, abs(HA) anggles greater
    than correct[0] are given a flat value correct[1].
    """

    dec = np.atleast_1d(dec)
    ha = np.atleast_1d(ha) % 360.

    if ha > 180:
        ha -= 360

    airmass = (
        np.sin(lat * np.pi / 180.) * np.sin(dec * np.pi / 180.) +
        np.cos(lat * np.pi / 180.) * np.cos(dec * np.pi / 180.) * np.cos(ha * np.pi / 180.))**(-1)

    if correct is not None:
        airmass[np.abs(ha) > correct[0]] = correct[1]

    if len(airmass) == 1:
        return airmass[0]
    else:
        return airmass


def mark_plate_complete(plate):
    """Sets the plugging status of ``plate`` to ``Complete``."""

    totoroDB = getConnection()
    plateDB = totoroDB.plateDB
    session = totoroDB.Session()

    plugging_status = [plugging.status.label for plugging in plate.pluggings]
    if 'Good' in plugging_status or len(plugging_status) == 0:
        return

    sorted_plugging = sorted(plate.pluggings, key=lambda plug: (plug.fscan_mjd, plug.fscan_id))

    last_plugging = sorted_plugging[-1]

    with session.begin():
        good_pk = session.query(
            plateDB.PluggingStatus.pk).filter(plateDB.PluggingStatus.label == 'Good').scalar()
        last_plugging.plugging_status_pk = good_pk

    return True


def isPlateComplete(plate,
                    format='plate_id',
                    forceCheckCompletion=False,
                    write_apocomplete=True,
                    overwrite=False,
                    mark_complete=True,
                    **kwargs):
    """Returns True if a plate is complete using the MaNGA logic.

    If ``forceCheckCompletion`` is False and the plugging is marked as
    complete, no plateCompletion check is performed (this saves some time). If
    ``write_apocomplete=True`` and the plate is complte, the apocomplete file
    will be written. If ``mark_complete=True``, changes the plugging status to
    ``Complete``.

    """

    from Totoro.dbclasses.plate import Plate

    if not isinstance(plate, Plate):
        if format.lower() not in ['pk', 'plate_id']:
            raise exceptions.TotoroError('format must be plate_id or pk.')
        plate = Plate(plate, format=format.lower(), **kwargs)

    plugStatus = [plugging.status.label for plugging in plate.pluggings]
    field_name = plate.field_name

    if 'Good' in plugStatus or 'Overridden Good' in plugStatus:
        plugComplete = True
    elif u'Overridden Incomplete' in plugStatus:
        plugComplete = False
    else:
        plugComplete = None

    if plugComplete is not None and forceCheckCompletion is False:
        if plugComplete is True and len(plate.getMockExposures()) == 0:
            if write_apocomplete:
                getAPOcomplete([plate], createFile=True, overwrite=overwrite)
            if mark_complete:
                mark_plate_complete(plate)
        return plugComplete

    completion_threshold = config['SN2thresholds']['completionThreshold']
    if plate.getPlateCompletion(includeIncompleteSets=False) >= completion_threshold:
        plateComplete = True
    else:
        plateComplete = False

    # Special plates can be completed in a number of sets, even if the
    # SN2 level are lower than the threshold.
    if not plateComplete:
        if field_name is not None and field_name in config['specialPrograms']:
            field_config = config['specialPrograms'][field_name]
            if 'complete_with_n_sets' in field_config:
                n_sets = len([ss for ss in plate.sets if ss.complete is True])
                if n_sets >= field_config['complete_with_n_sets']:
                    plateComplete = True

    if (plateComplete is True and np.isnan(plate.getCumulatedSN2(includeIncomplete=False)).any()):
        log.debug('plate_id={0}: not all cameras have been correctly reduced. '
                  'Setting plateComplete=False.'.format(plate.plate_id))
        plateComplete = False

    if plugComplete is not None:
        if plugComplete is not plateComplete:
            warnings.warn(
                'plate={0}: plugging status is {1} but calculated '
                'status is {2}.'.format(plate.plate_id, 'complete'
                                        if plugComplete else 'incomplete', 'complete'
                                        if plateComplete else 'incomplete'),
                exceptions.TotoroUserWarning)

    completion_status = plugComplete or plateComplete

    if completion_status is True and len(plate.getMockExposures()) == 0:
        if write_apocomplete:
            getAPOcomplete([plate], createFile=True, overwrite=overwrite)
        if mark_complete:
            mark_plate_complete(plate)

    return completion_status


def getAPOcomplete(plates,
                   format='plate_id',
                   SN2_blue=None,
                   SN2_red=None,
                   limitSN=False,
                   func=np.max,
                   createFile=False,
                   reindex_sets=True,
                   **kwargs):
    """Returns a dictionary with the APOcomplete output.

    Parameters
    ----------
    plates : list of `Totoro.Plate` instances or list of ints
        Either a single `Totoro.Plate` instance or and integer, or a list of
        them. If an integer (or list), the appropriate plate(s) will be
        obtained from that value and the `format` paramenter.
    format : string
        If `plates` are integers, the field of the plate table on which to
        perform the query. Normally either `'plate_id'` or `'pk'`.
    SN2_blue, SN2_red : None or float
        The SN2 plate thresholds in blue and red, respectively. If None,
        the values in config.SN2thresholds.plate(Blue|Red) will be used.
    limitSN : bool
        If True, the function will use only the combination of sets that gets
        a cumulated SN2 closer (but higher) than `SN2_blue` and `SN2_red.
        If the sum of all the valid sets is not enough to reach the SN2
        thresholds, all the sets will be used. That is, if a plate has four
        valid sets but three of them are enough to reach the SN2 thresholds,
        only those are used.
    func : function
        If limitSN is True and several combinations of sets meet the
        requirement of having SN2 higher than the thresholds, this function
        is used to determine which combination to use. Usual options are
        `np.max` to use the combination that gives higher SN2 with fewer sets,
        or `np.min` to use the combination that gets closer to the SN2
        thresholds (but still above them).
    createFile : bool
        If True, `createAPOcompleteFile` is called for each of the plates.
    reindex_sets : bool
        If ``True``, reindexes the ``set_pk`` to ``1, 2, ...``.
    kwargs : dict
        Additional parameters to be passed to `Totoro.Plates` and to
        `createAPOcompleteFile`.

    """

    from Totoro.dbclasses import Plate

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

        if isPlateComplete(plate, write_apocomplete=False) is False:
            warnings.warn(
                'plate_id={0} is not complete. APOcomplete output '
                'must not be used.'.format(plate.plate_id), exceptions.TotoroUserWarning)

        APOcomplete[plate.plate_id] = []

        validSets = plate.getValidSets()

        if limitSN:
            for nSets in range(1, len(validSets) + 1):

                combSets = list(combinations(validSets, nSets))
                SN2 = np.array([_cumulatedSN2(sets) for sets in combSets])

                overSN2 = [(combSets[ii], SN2[ii]) for ii in range(len(combSets))
                           if SN2[ii][0] >= SN2_blue and SN2[ii][1] >= SN2_red]

                if len(overSN2) > 0:

                    relativeSN2 = np.array(
                        [(overSN2[ii][1][0] / SN2_blue) * (overSN2[ii][1][1] / SN2_red)
                         for ii in range(len(overSN2))])

                    setsToAPOcomplete = overSN2[np.where(
                        relativeSN2 == func(relativeSN2))[0][0]][0]

                    break

        if setsToAPOcomplete is None:
            setsToAPOcomplete = validSets

        for set_ii, ss in enumerate(setsToAPOcomplete):
            for exp in ss.totoroExposures:

                mjd = exp.getMJD()
                pk = ss.pk if reindex_sets is False else (set_ii + 1)
                dPos = exp.ditherPosition.upper()
                nExp = exp.exposure_no

                APOcomplete[plate.plate_id].append([plate.plate_id, mjd, pk, dPos, nExp])

        if len(setsToAPOcomplete) > 0:
            apoCompleteSN2 = _cumulatedSN2(setsToAPOcomplete)
        else:
            apoCompleteSN2 = np.array([0., 0.])

        if apoCompleteSN2[0] >= SN2_blue and apoCompleteSN2[1] >= SN2_red:
            log.info('APOcomplete for plate_id={0} returned with '
                     'SN2_blue={1:.1f}, SN2_red={2:.1f}.'.format(plate.plate_id, apoCompleteSN2[0],
                                                                 apoCompleteSN2[1]))
        else:
            warnings.warn(
                'plate_id={0} has SN2_blue={1:.1f}, SN2_red={2:.1f},'
                ' which is lower than the thresholds.'.format(plate.plate_id, apoCompleteSN2[0],
                                                              apoCompleteSN2[1]),
                exceptions.TotoroUserWarning)

    if createFile:
        createAPOcompleteFile(APOcomplete, **kwargs)

    return APOcomplete


def createAPOcompleteFile(APOcomplete, path=None, overwrite=False, svn_add=True):
    """Writes the APOcomplete file in Yanny format."""

    for plate in APOcomplete:

        plateXX = '{:06d}'.format(plate)[0:4] + 'XX'
        default_path = os.path.join(os.environ['MANGACORE_DIR'], 'apocomplete', plateXX)

        path = default_path if path is None else path

        apocompPath = os.path.join(path, 'apocomp-{0:04d}.par'.format(plate))

        if os.path.exists(apocompPath):
            if overwrite:
                warnings.warn('apocomplete path {} exists but '
                              'overwriting it.'.format(path), exceptions.TotoroUserWarning)
            else:
                log.debug('apocomplete path {} exists; not ' 'overwriting it.'.format(path))
                return

        if not os.path.exists(path):
            os.makedirs(path)

        data = APOcomplete[plate]

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
            expstr = '{0} {1} {2} {3} {4} {5}\n'.format('APOCOMP', plateid[ii], mjd[ii], setno[ii],
                                                        dpos[ii], expno[ii])
            ff.write(expstr)

        ff.close()

        if svn_add:
            try:
                os.chdir(path)
                result = subprocess.call('svn add {}'.format(apocompPath), shell=True)
                if result > 0:
                    warnings.warn('svn add {} failed with error {}'.format(apocompPath, result),
                                  exceptions.TotoroUserWarning)
                    return
            except Exception:
                warnings.warn('svn add {} failed with unknown error'.format(apocompPath),
                              exceptions.TotoroUserWarning)
                return

    return path


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

    totoroDB = getConnection()
    plateDB = totoroDB.plateDB
    session = totoroDB.session

    from Totoro.dbclasses import Plate

    if isinstance(plate, (Plate, plateDB.Plate)):
        pass
    else:
        try:
            with session.begin():
                plate = session.query(plateDB.Plate).filter(plateDB.Plate.plate_id == plate).one()
        except Exception:
            return False

    for survey in plate.surveys:
        if (survey.label == 'MaNGA' and plate.currentSurveyMode is not None and
                plate.currentSurveyMode.label in ['MaNGA dither', 'MaNGA 10min']):
            return True

    return False


def checkOpenSession():
    """Raises an error if Totoro is being run from inside an open session."""

    totoroDB = getConnection()
    session = totoroDB.Session()

    try:
        with session.begin():
            session.commit()
    except ResourceClosedError:
        pass
    except InvalidRequestError as ee:
        if 'A transaction is already begun' in str(ee):
            raise exceptions.TotoroSubtransactionError(
                'Totoro is being run within an open SQLalchemy session. '
                'Please, modify your code to avoid this.')
        else:
            raise exceptions.TotoroSubtransactionError(
                'Failed while checking session status. Error message is: {0}'.format(str(ee)))
    except Exception as ee:
        raise exceptions.TotoroSubtransactionError(
            'Failed while checking session status. Error message is: {0}'.format(str(ee)))


def get_closest_holes(plateid):
    """Calculates the minimum distance between holes in the plateHoles.

    Returns a tuple containing the distance between the closest pair
    (in mm), the xFocal and yFocal of the two closest holes, and their
    hole types.

    It checks MANGA, OBJECT, and GUIDE holes.
    We ignore the separation between two OBJECT holes since those
    correspond always APOGEE or eBOSS fibres.

    """

    # The holeType to take into account when finding the closest pair`
    valid_holes = ['GUIDE', 'MANGA', 'OBJECT']

    if 'PLATELIST_DIR' not in os.environ:
        raise ValueError('cannot access the platelist product')

    plate6 = '{0:06d}'.format(plateid)
    short_platedir = plate6[0:4] + 'XX'
    plateHoles_path = os.path.join(os.environ['PLATELIST_DIR'], 'plates', short_platedir, plate6,
                                   'plateHoles-{0}.par'.format(plate6))

    if not os.path.exists(plateHoles_path):
        plateHoles_path = plateHoles_path.replace('plateHoles', 'plateHolesSorted')
        if not os.path.exists(plateHoles_path):
            raise ValueError('cannot find plateHoles for plate {0}'.format(plateid))

    plateHoles = yanny(plateHoles_path)['STRUCT1']

    mask = np.in1d(plateHoles['holetype'], valid_holes)
    holes = plateHoles[mask]

    focal = np.zeros((len(holes['xfocal']), 2), dtype=np.float64)
    focal[:, 0] = holes['xfocal']
    focal[:, 1] = holes['yfocal']

    distances = pdist(focal)
    pdist_indices = list(combinations(list(range(focal.shape[0])), 2))

    for kk in np.argsort(distances):
        ii, jj = pdist_indices[kk]
        hole1_type = holes['holetype'][ii]
        hole2_type = holes['holetype'][jj]
        if hole1_type != 'OBJECT' or hole2_type != 'OBJECT':
            return (distances[kk], focal[ii], focal[jj], hole1_type, hole2_type)


def avoid_cart_2(plate):
    """Finds closest pair of holes and decides whether to use cart 2.

    This is because cart 2 has heat shrinks around some of the fibres,
    which effectively increases their ferrule size. If the distance between
    holes is smaller than the enlarged ferrule sizes, we should avoid
    cart 2, if possible.

    The results are cached to avoid having to reopen the yanny file.

    """

    if plate in _avoid_cart2_cahche:
        return _avoid_cart2_cahche[plate]

    distance, hole1_focal, hole2_focal, \
        hole1_type, hole2_type = get_closest_holes(plate)

    min_distance = 0.5 * (config['ferruleSizes'][hole1_type] + config['ferruleSizes'][hole2_type])

    if distance <= min_distance:
        _avoid_cart2_cahche[plate] = True
        return True
    else:
        _avoid_cart2_cahche[plate] = False
        return False
