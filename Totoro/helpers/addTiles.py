# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides some functions and classes to add tiles to
mangaDB from a tiling catalogue.
"""

from __future__ import print_function
from ..plateDB.dataModel import MangaDB_Field, db, session
# from ..exceptions import TotoroError
from astropy import table
import numpy as np

engine = db.engine

__all__ = ['addFromTilingCatalogue']

# Grouping values for different indexes and keys
TILE_ID_ROOT = 0


def addFromTilingCatalogue(catalogue, isSample=True):
    """Reads a tiling catalogue and adds the information to plateDB.

    This function reads a sample or tiling catalogue and uptates mangaDB.field.

    """

    if not isSample:
        catData = readTilingCatalogue(catalogue)
    else:
        catData = readSampleCatalogue(catalogue)

    session.rollback()

    for tile in catData:
        addTiles(catData, commit=False)
    session.commit()


def addTiles(tiles, commit=True):
    """Adds fields to the database."""

    data = []

    for tile in tiles:

        tileID = TILE_ID_ROOT + tile['LOCATIONID']
        RA = tile['PLATERA']
        Dec = tile['PLATEDEC']
        expectedNOfVisits = 1
        shared = False
        name = 'Tile {0:d}'.format(int(tileID))

        dd = dict(
            name=name, center_ra=RA, center_dec=Dec,
            location_id=tileID,
            expected_no_visits=expectedNOfVisits,
            shared=shared)

        data.append(dd)

    engine.execute(
        MangaDB_Field.__table__.insert(), data)

    if commit:
        session.commit()

    return


def readTilingCatalogue(catalogue, id='ID', ra='RA', dec='DEC'):
    """Reads a tiling catalogue in FITS format."""

    cat = table.Table.read(catalogue, format='fits')
    if id != 'ID':
        cat.rename_column(id, 'ID')
    if ra != 'RA':
        cat.rename_column(ra, 'RA')
    if dec != 'Dec':
        cat.rename_column(dec, 'Dec')
    return cat


def readSampleCatalogue(catalogue):

    sample = table.Table.read(catalogue, format='fits')
    locationIDs, idx = np.unique(sample['LOCATIONID'],
                                 return_index=True)
    validIdx = locationIDs > 0
    locationID_Idx = idx[validIdx]

    tabPlates = sample[locationID_Idx]['LOCATIONID', 'PLATERA', 'PLATEDEC']
    tabPlates.sort('LOCATIONID')

    return tabPlates
