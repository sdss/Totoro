#!/usr/bin/env python
# encoding: utf-8
"""
loadMangaPlates.py

Created by José Sánchez-Gallego on 6 Mar 2015.
Licensed under a 3-clause BSD license.

Revision history:
    6 Mar 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from Totoro import config, readPath
from Totoro.db import getConnection
from Totoro.exceptions import TotoroError
from pydl.pydlutils.yanny import yanny
from astropy import table
import glob
import os
import sys

db = getConnection()
session = db.Session()


def getMangaTileIDs():
    """Returns a dictionary of {plateid: manga_tileid} from the
    plateTargets files stored in mangacore."""

    mangacorePath = readPath(config['fields']['mangacore'])
    plateTargets = glob.glob(
        os.path.join(mangacorePath,
                     'platedesign/platetargets/plateTargets-*.par'))

    if len(plateTargets) == 0:
        raise TotoroError('no plateTargets files found.')

    mangaTileIDs = {}
    neverobserve = {}
    for plateTargetsFile in plateTargets:
        pT = yanny(plateTargetsFile)['PLTTRGT']
        for target in pT:
            if (target['plateid'] not in mangaTileIDs and
                    target['manga_tileid'] > 0):
                mangaTileIDs[target['plateid']] = target['manga_tileid']
            if target['plateid'] not in neverobserve:
                neverobserve[target['plateid']] = target['neverobserve']

    return mangaTileIDs, neverobserve


def readSpecialPlates():
    """Returns an astropy.table.Table instance with the data for special
    plates."""

    specialPlatesFile = os.path.join(os.path.dirname(__file__),
                                     '../data/specialPlates.dat')

    return table.Table.read(specialPlatesFile, format='ascii.commented_header',
                            delimiter=';')


def loadMangaPlates():
    """Loads plate information into mangaDB.Plate"""

    mangaTileIDs, neverobserve = getMangaTileIDs()

    specialPlates = readSpecialPlates()
    specialPlates['comment'].fill_value = ''
    specialPlates = specialPlates.filled()

    with session.begin():
        allPlates = session.query(db.plateDB.Plate).join(
            db.plateDB.PlateToSurvey, db.plateDB.Survey,
            db.plateDB.SurveyMode).filter(
                db.plateDB.Survey.label == 'MaNGA',
                db.plateDB.SurveyMode.label.like('%MaNGA%')).all()

    with session.begin():

        for nn, plate in enumerate(allPlates):

            try:
                newPlate = session.query(db.mangaDB.Plate).filter(
                    db.mangaDB.Plate.platedb_plate_pk == plate.pk).one()
            except:
                newPlate = db.mangaDB.Plate()

            if plate.plate_id in mangaTileIDs:
                newPlate.manga_tileid = mangaTileIDs[plate.plate_id]
            else:
                newPlate.manga_tileid = None

            if plate.plate_id in specialPlates['plateid']:
                row = specialPlates[specialPlates['plateid'] ==
                                    plate.plate_id]
                newPlate.special_plate = True
                newPlate.all_sky_plate = bool(row['all_sky_plate'][0])
                newPlate.commissioning_plate = bool(
                    row['commissioning_plate'][0])
                newPlate.comment = row['comment'][0]
            else:
                newPlate.special_plate = False
                newPlate.all_sky_plate = False
                newPlate.commissioning_plate = False
                newPlate.comment = ''

            if plate.plate_id in neverobserve:
                newPlate.neverobserve = bool(neverobserve[plate.plate_id])
            # Some harcoded values that may not appear in neverobserve
            elif plate.plate_id in [7566, 7567, 7568]:
                newPlate.neverobserve = True
            else:
                newPlate.neverobserve = False

            newPlate.platedb_plate_pk = plate.pk

            session.add(newPlate)

            sys.stdout.write('\rLoading plates: {0:.0f}%'
                             .format(nn / float(len(allPlates)) * 100.))
            sys.stdout.flush()
    return


if __name__ == '__main__':
    loadMangaPlates()
