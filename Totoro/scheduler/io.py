#!/usr/bin/env python
# encoding: utf-8
"""

io.py

Created by José Sánchez-Gallego on 19 Oct 2015.
Licensed under a 3-clause BSD license.

Revision history:
    19 Oct 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division, print_function

import os

import numpy as np
from astropy import table

from Totoro.dbclasses import Exposure, Field, Set


def saveExposures(plates, outfile, startDate=None):
    """Writes information about simulated plates to a FITS file."""

    template = table.Table(
        None,
        names=[
            'plate_id', 'manga_tileid', 'set_pk', 'real_set', 'start_jd', 'exposure_time',
            'dither_position', 'ra', 'dec', 'sn2values'
        ],
        dtype=[int, int, int, int, float, float, 'S1', float, float, (float, 4)])

    for plate in plates:
        setID = 1
        for ss in plate.sets:
            if not ss.isMock:
                set_pk = ss.pk
                real_set = 1
            else:
                set_pk = setID
                real_set = 0
                setID += 1

            if plate.isMock:
                plate_id = -999
            else:
                plate_id = plate.plate_id

            ra = plate.ra
            dec = plate.dec

            manga_tileid = (plate.manga_tileid if plate.manga_tileid is not None else -999)

            for exp in ss.totoroExposures:

                if not exp.isMock:
                    continue

                template.add_row((plate_id, manga_tileid, set_pk, real_set,
                                  exp.getJD()[0], exp.exposure_time, exp.ditherPosition, ra, dec,
                                  exp.getSN2Array()))

    # Records the start date
    if startDate:
        template.meta['STRTDATE'] = startDate
    else:
        template.meta['STRTDATE'] = np.min(template['start_jd'])

    if os.path.exists(outfile):
        os.remove(outfile)

    template.write(outfile, format='fits')


def createExposure(row):
    """Returns a mock exposure."""

    exp = Exposure.createMockExposure(
        startTime=row['start_jd'],
        expTime=row['exposure_time'],
        ditherPosition=row['dither_position'],
        ra=row['ra'],
        dec=row['dec'],
        sn2values=row['sn2values'],
        seeing=1.0)

    return exp


def removeExposures(plates, startDate):
    """Removes all real exposures taking after startDate."""

    for plate in plates:
        for ss in plate.sets:
            ss.isMock = True
            validExps = [exp for exp in ss.totoroExposures if exp.getJD()[0] < startDate]
            if len(validExps) != len(ss.totoroExposures):
                plate._modified = True
                plate._useOnlyCompletion = True
            ss.totoroExposures = validExps
        plate.sets = [ss for ss in plate.sets if len(ss.totoroExposures) > 0]

    return plates


def restoreExposures(exposureFile, plates=[]):
    """Restores exposures to a list of plates."""

    data = table.Table.read(exposureFile)

    if 'STRTDATE' in data.meta:
        startDate = data.meta['STRTDATE']
        if startDate is not None and startDate > 0:
            plates = removeExposures(plates, startDate)

    for plate in plates:
        plate_id = plate.plate_id
        exps = data[data['plate_id'] == plate_id]

        realSetExps = exps[exps['real_set'] == 1]
        for exp in realSetExps:
            mockExp = createExposure(exp)
            for ss in plate.sets:
                if ss.pk == exp['set_pk']:
                    ss.totoroExposures.append(mockExp)

        mockSetExps = exps[exps['real_set'] == 0]
        uniqueSets = np.unique(mockSetExps['set_pk'])
        for ss in uniqueSets:
            setExps = [createExposure(exp) for exp in mockSetExps[mockSetExps['set_pk'] == ss]]
            newSet = Set.fromExposures(setExps)
            plate.sets.append(newSet)

    fields = data[data['plate_id'] == -999.]
    uniqueMangaTileIDs = np.unique(fields['manga_tileid'])
    for mangaTileID in uniqueMangaTileIDs:
        exps = fields[fields['manga_tileid'] == mangaTileID]
        newPlate = Field.createMockPlate(ra=exps[0]['ra'], dec=exps[0]['dec'])
        newPlate.manga_tileid = mangaTileID
        uniqueSets = np.unique(exps['set_pk'])
        for ss in uniqueSets:
            setExps = [createExposure(exp) for exp in exps[exps['set_pk'] == ss]]
            newSet = Set.fromExposures(setExps)
            newPlate.sets.append(newSet)
        plates.append(newPlate)

    return plates
