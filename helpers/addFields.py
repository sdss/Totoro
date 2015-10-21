#!/usr/bin/env python
# encoding: utf-8
"""
addFields.py

Created by José Sánchez-Gallego on 16 Apr 2014.
Licensed under a 3-clause BSD license.

A quick routine to add the fields from a tiling catalogue. This is for
testing purposes and the final field management must be done using the
designated product to manage the MaNGA sample.

Revision history:
    16 Apr 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from ..APOplateDB import sampleMangaDB, Session, engine
from astropy import table


session = Session()


def addFields(catalogue):
    """Accepts a catalogue in FITS format and updates the Field table.

    Note that the Field table is purged before proceeding.

    """

    with session.begin(subtransactions=True):
        session.query(sampleMangaDB.Field).delete()

    tiles = table.Table.read(catalogue)

    fieldData = []

    for tile in tiles:
        fieldData.append(
            {'location_id': tile['ID'],
             'center_ra': tile['RA'],
             'center_dec': tile['DEC'],
             'expected_no_visits': 1})

    with session.begin(subtransactions=True):
        engine.execute(sampleMangaDB.Field.__table__.insert(), fieldData)
