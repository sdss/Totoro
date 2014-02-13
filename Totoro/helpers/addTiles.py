# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides some functions and classes to add tiles to
mangaDB from a tiling catalogue.
"""

from __future__ import print_function
from ..plateDB.dataModel import *
from ..exceptions import TotoroError


__all__ = ['addFromTilingCatalogue']


# Grouping values for different indexes and keys
TILE_ID_ROOT = 20000


def addFromTilingCatalogue(catalogue):
    """Reads a tiling catalogue and adds the information to plateDB.

    This function reads a tiling catalogue and uptates plateDB to include
    the tiles as plates in design phase.

    """

    session.rollback()
    catData = readTilingCatalogue(catalogue)

    for tile in catData:
        addTile(tile, commit=False)
    session.commit()


def addTile(tile, commit=True):
    """Adds a single tile to the database."""

    tileID = TILE_ID_ROOT + tile['ID']
    RA = tile['RA']
    Dec = tile['Dec']

    qq = session.query(MangaDB_Tile).filter(MangaDB_Tile.id == tileID)

    if qq.count() == 1:
        qq.update(dict(id=tileID, ra_centre=RA, dec_centre=Dec))

    elif qq.count() == 0:
        session.add(MangaDB_Tile(id=tileID, ra_centre=RA, dec_centre=Dec))
    else:
        raise TotoroError('Trying to add tile ID {0} but '.format(tileID) +
                          'more than one row in mangaDB.tile ' +
                          'already have that ID.')

    if commit:
        session.commit()

    return


def readTilingCatalogue(catalogue, id='ID', ra='RA', dec='DEC'):
    """Reads a tiling catalogue in FITS format."""
    from astropy import table
    cat = table.Table.read(catalogue, format='fits')
    if id != 'ID':
        cat.rename_column(id, 'ID')
    if ra != 'RA':
        cat.rename_column(ra, 'RA')
    if dec != 'Dec':
        cat.rename_column(dec, 'Dec')
    return cat
