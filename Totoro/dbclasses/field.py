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
import Totoro
import plate as TotoroPlate
from Totoro.db import getConnectionFull
from Totoro.exceptions import TotoroError
import os
from astropy import table
import numpy as np
from numbers import Integral


__all__ = ['Fields', 'Field', 'getTilingCatalogue']


def getTilingCatalogue(tilingCatalogue=None):
    """Returns the tiling catalogue"""

    tilingCatalogue = Totoro.readPath(
        Totoro.config['fields']['tilingCatalogue']) \
        if tilingCatalogue is None else Totoro.readPath(
            tilingCatalogue)

    if not os.path.exists(tilingCatalogue):
            raise TotoroError('tiling catalogue {0} does not exist.'
                              .format(os.path.realpath(tilingCatalogue)))

    tiles = table.Table.read(tilingCatalogue)

    return tiles


class Fields(list):

    def __init__(self, tilingCatalogue=None, rejectDrilled=True, silent=False,
                 **kwargs):
        """Returns a list of `Totoro.Field` instances.

        Reads a tiling catalogue and returns a list of `Totoro.Field` objects,
        one for each tile. If `rejectDrilled=True`, tiles that have already
        been drilled are skipped.

        Parameters
        ----------
        tilingCatalogue : str or None
            The path to the tiling catalogue to be used
        rejectDrilled : bool
            If True, already drilled tiles are skipped.
        silent : bool
            If True, does limited logging.
        kwarg : dict
            Additional arguments to be passed during `Field` creation.

        """

        self._tiles = getTilingCatalogue(tilingCatalogue=tilingCatalogue)

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

    def _rejectDrilled(self, silent=False):
        """Rejects plates in self that have been drilled."""

        __, Session, plateDB, __ = getConnectionFull()
        session = Session()

        with session.begin():
            plates = session.query(plateDB.Plate).join(
                plateDB.PlateToSurvey, plateDB.Survey, plateDB.SurveyMode
            ).filter(plateDB.Survey.label == 'MaNGA',
                     plateDB.SurveyMode.label == 'MaNGA dither').order_by(
                         plateDB.Plate.plate_id).all()

        plateMangaTileIds = np.unique(
            [plate.mangadbPlate.manga_tileid
             for plate in plates if plate is not None and
             plate.mangadbPlate is not None])

        drilledFields = []
        for field in self:
            if field.manga_tileid in plateMangaTileIds:
                drilledFields.append(field)

        nRemoved = len(drilledFields)
        for dField in drilledFields:
            self.remove(dField)

        logMsg = ('rejected {0} fields because they have already been drilled'
                  .format(nRemoved))
        if not silent:
            Totoro.log.info(logMsg)
        else:
            Totoro.log.debug(logMsg)

    def removeField(self, inp):
        """Removes a field."""

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


class Field(TotoroPlate.Plate):

    def __repr__(self):
        return '<Field: manga_tileid={0}>'.format(self.manga_tileid)

    @property
    def isComplete(self):
        if self._complete is not None:
            return self._complete
        else:
            if self.getPlateCompletion(includeIncompleteSets=False) >= 1.:
                self._complete = True
                return True
            else:
                return False
