#!/usr/bin/env python
# encoding: utf-8
"""
field.py

Created by José Sánchez-Gallego on 17 Feb 2014.
Licensed under a 3-clause BSD license.

Revision history:
    17 Feb 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.internal.manga import Totoro
import plate
from sdss.internal.manga.Totoro.exceptions import TotoroError
import os
from astropy import table
import numpy as np
from numbers import Integral


db = Totoro.TotoroDBConnection()
session = db.Session()


class Fields(list):

    def __init__(self, tilingCatalogue=None, rejectDrilled=True, silent=False,
                 **kwargs):

        self.tilingCatalogue = Totoro.readPath(
            Totoro.config['fields']['tilingCatalogue']) \
            if tilingCatalogue is None else Totoro.readPath(
                tilingCatalogue)

        self._tiles = self._getTiles(**kwargs)

        mockFields = []
        for tile in self._tiles:
            mockField = Field.createMockPlate(
                ra=tile['RA'], dec=tile['DEC'], silent=True, **kwargs)
            mockField.manga_tileid = int(tile['ID'])
            mockFields.append(mockField)

        list.__init__(self, mockFields)

        logMsg = 'loaded {0} fields from tiling catalogue'.format(len(self))
        if not silent:
            Totoro.log.info(logMsg)
        else:
            Totoro.log.debug(logMsg)

        if rejectDrilled:
            self._rejectDrilled(silent=silent, **kwargs)

    def _getTiles(self, **kwargs):

        if not os.path.exists(self.tilingCatalogue):
            raise TotoroError('tiling catalogue {0} does not exist.'
                              .format(os.path.realpath(self.tilingCatalogue)))

        tiles = table.Table.read(self.tilingCatalogue)

        return tiles

    def _rejectDrilled(self, silent=False, rejectMode='manga_tileid',
                       **kwargs):
        """Rejects plates in self that have been drilled. rejectMode is a
        placeholder for future functionality when fields might be rejected
        based on the coordinates of drilled plates."""

        allDrilled = plate.Plates.getAll(onlyIncomplete=False, silent=True,
                                         updateSets=False)
        mangaTileIDs = map(lambda xx: xx.getMangaTileID(), allDrilled)

        nRemoved = 0
        for mangaTileID in mangaTileIDs:
            result = self.removeField(mangaTileID)
            if result:
                nRemoved += 1

        logMsg = ('rejected {0} fields because they have already been drilled'
                  .format(nRemoved))
        if not silent:
            Totoro.log.info(logMsg)
        else:
            Totoro.log.debug(logMsg)

    def removeField(self, inp):

        if isinstance(inp, Integral):
            mangaTileIDs = np.array([ff.manga_tileid for ff in self], int)
            idx = np.where(mangaTileIDs == inp)[0]

            if idx.size != 1:
                return None
            else:
                self.remove(self[idx[0]])
                return True

        elif isinstance(inp, Field):
            self.remove(inp)
            return True

        elif inp is None:
            return None

        else:
            raise TotoroError('input must be a manga_tileid integer or a Field'
                              ' instance.')


class Field(plate.Plate):

    def __repr__(self):
        return '<Field: manga_tileid={0:d}>'.format(self.manga_tileid)
