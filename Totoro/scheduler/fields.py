#!/usr/bin/env python
# encoding: utf-8
"""
fields.py

Created by José Sánchez-Gallego on 17 Feb 2014.
Licensed under a 3-clause BSD license.

Revision history:
    17 Feb 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from ..core.defaults import SAMPLE, LATITUDE, EXPTIME, NDITHERS
from ..core.defaults import AVG_SN_RED, AVG_SN_BLUE, ALPHA_RED, ALPHA_BLUE, \
    R_SN2, B_SN2, EFFICIENCY, SKY_PRIORITY, COMPLETION_PRIORITY
from ..exceptions import TotoroError, TotoroWarning
import warnings
import os
import numpy as np
from astropy import table as tt
from ..utils import getCompletedFields
from ..plateDB.dataModel import session, MangaDB_Field
from astropy import coordinates as coo
from astropy import units as uu
from ..utils import DustMap
from .exposure import Exposure


try:
    DUST_MAP = DustMap()
except:
    DUST_MAP = None


class Fields(list):

    def __init__(self, table=False, sample=False, **kwargs):

        if table is not False:
            self._fields = table
            self.dataSource = 'table'
        elif sample is not False:
            self._fields = self._fromSample(sample)
            self.dataSource = 'sample'
        else:
            self._fields = self._fromMangaDB()
            self.dataSource = 'mangaDB'

        list.__init__(
            self,
            [Field(field) for field in self._fields]
        )

    def _fromSample(self, sample):
        if isinstance(sample, basestring):
            pass
        else:
            if SAMPLE() is not None:
                sample = SAMPLE()
            elif 'MANGASAMPLE' in os.environ:
                sample = os.environ['MANGASAMPLE']
            else:
                raise TotoroError('no sample source.')

        return self._getFieldsFromSample(sample)

    @staticmethod
    def _getFieldsFromSample(sample):
        """Returns a table with the fields in a sample catalogue.

        This method returns a `Table` with LOCATIONID, PLATERA and PLATEDEC
        for each unique locationID in a sample catalogue.

        """

        if isinstance(sample, basestring):
            sampleData = tt.Table.read(sample)
        elif isinstance(sample, tt.Table):
            sampleData = sample
        else:
            raise TotoroError('sample format not understood.')

        locationIDs, idx = np.unique(sampleData['LOCATIONID'],
                                     return_index=True)
        validIdx = locationIDs > 0
        locationID_Idx = idx[validIdx]
        return sampleData[locationID_Idx]['LOCATIONID', 'PLATERA', 'PLATEDEC']

    @staticmethod
    def _fromMangaDB():

        fields = session.query(
            MangaDB_Field.location_id, MangaDB_Field.center_ra,
            MangaDB_Field.center_dec, MangaDB_Field.pk)

        cols = zip(*fields.all())

        table = tt.Table(
            cols, names=['LOCATIONID', 'PLATERA', 'PLATEDEC', 'FIELD_PK'],
            dtype=[int, float, float, int])

        return table

    def removeCompleted(self):

        if self.dataSource == 'table' or \
                self.dataSource == 'sample':
            warnings.warn('Totoro.Tiles.removeCompleted inly works '
                          'when dataSource=\'mangaDB\'. No fields have '
                          'been removed.', TotoroWarning)
        else:
            completedFields = getCompletedFields()
            for cc in completedFields:
                self.deleteField(cc)

    def deleteField(self, locationID):
        locIDs = [ff.locationID for ff in self]
        if locationID not in locIDs:
            raise TotoroError('locationID not found.')
        idx = locIDs.index(locationID)
        self.pop(idx)

    def getPriorities(self, sortedIndices=False):
        """Returns a list with the priorities of the fields.

        This method returns a list with the priority of each field in the
        object. If ``sortIndices=True``, a list of indices sorted from higher
        to lower priority is also returned.

        """

        priorities = [ff.getPriority() for ff in self]
        order = np.argsort(priorities)[::-1]

        if sortedIndices:
            return (priorities, order)
        else:
            return priorities


class Field(object):

    def __init__(self, field):

        self._priority = None

        self.locationID = int(field['LOCATIONID'])
        self.plateCentre = coo.ICRS(ra=field['PLATERA'], dec=field['PLATEDEC'],
                                    unit=(uu.degree, uu.degree))

        if 'FIELD_PK' in field.colnames:
            self.fieldPK = field['FIELD_PK']

        self.HAlimits = self.mlhalimit()

        if DUST_MAP is None:
            self.iIncrease = 1
            self.gIncrease = 1
        else:
            self.gIncrease, self.iIncrease = DUST_MAP.eval(
                self.plateCentre.ra.deg, self.plateCentre.dec.deg)

        self.SN2_Blue = None
        self.SN2_Red = None

        self.nExposuresH0 = self.getCompletionTime(HA=0.)

        self.exposures = []
        self.isCompleted = False
        self.setCompleted = True

    @property
    def nExposures(self):
        return len(self.exposures)

    @property
    def observabilityWindow(self):
        if self.setCompleted is True:
            return [-self.HAlimits, self.HAlimits]
        else:
            return [-self.HAlimits, self.HAlimits]  # To be changed

    def getPriority(self, LST):

        HA = LST.hour - self.plateCentre.ra.hour
        HA.wrap_at(24 * uu.hour, inplace=True)

        airMass = self.computeAirmass(HA.deg)
        if airMass > 2:
            return 0

        skyPriority = 1. / (airMass ** ALPHA_RED() * self.iIncrease)
        skyPriority *= SKY_PRIORITY()

        completionPriority = self.exposuresCompleted / self.nExposuresH0
        completionPriority *= COMPLETION_PRIORITY()

        return skyPriority + completionPriority

    def getCompletionTime(self, HA=0, direction='both'):

        haArray = [0]

        if direction == 'both':
            for ii in range(1, 10):
                haArray += [ii, -ii]
        elif direction == 'forward':
            for ii in range(1, 20):
                haArray += [ii]
        elif direction == 'backwards':
            for ii in range(1, 20):
                haArray += [-ii]

        haArray = np.array(haArray, dtype=float)
        haArray = HA + 15. * haArray * EXPTIME() / 60. / EFFICIENCY()

        airmasses = self.computeAirmass(haArray)

        snRed = self.SN2_Red + np.cumsum(
            AVG_SN_RED() / airmasses ** ALPHA_RED() / self.iIncrease)
        snBlue = self.SN2_Blue + np.cumsum(
            AVG_SN_BLUE() / airmasses ** ALPHA_BLUE() / self.iIncrease)

        if snBlue[-1] < B_SN2() or snRed[-1] < R_SN2():
            return None

        idxRed = np.where(snRed >= R_SN2())[0][0]
        idxBlue = np.where(snBlue >= B_SN2())[0][0]
        idx = np.max([idxRed, idxBlue])

        nExp = idx + 1
        nExp = np.int(np.ceil(np.float(nExp) / NDITHERS()) * NDITHERS())

        return nExp

    def computeAirmass(self, ha, correct=[75., 10.]):

        lat = LATITUDE()
        dec = self.plateCentre.dec.deg

        airmass = (np.sin(lat * np.pi / 180.) * np.sin(dec * np.pi / 180.) +
                   np.cos(lat * np.pi / 180.) * np.cos(dec * np.pi / 180.) *
                   np.cos(ha * np.pi / 180.)) ** (-1)

        if correct is not None:
            if hasattr(airmass, '__getitem__'):
                airmass[np.abs(ha) > correct[0]] = correct[1]
            else:
                if np.abs(ha) > correct[0]:
                    airmass = correct[1]

        return airmass

    def mlhalimit(self):
        """Returns HA limits.

        Calculates the maximum HAs acceptable for a list of declinations.
        Uses the polinomial fit by David Law and a omega limit of 0.5.

        """

        funcFit = np.array([1.59349, 0.109658, -0.00607871,
                            0.000185393, -2.54646e-06, 1.16686e-08])[::-1]

        dec = self.plateCentre.dec.degree
        if dec < -10 or dec > 80:
            return 0.
        else:
            return np.abs(np.polyval(funcFit, dec))

    def addExposure(self, **kwargs):
        """Adds an exposure to the field."""
        self.exposures.append(Exposure(**kwargs))
        self._uptdateCompletion()

    # def _uptdateCompletion(self):

    def __repr__(self):
        return(
            '<Totoro.Field (designID={0:d}, RA={1:.4f}, Dec={2:.4f})>'.format(
                self.locationID, self.RA, self.Dec)
            )
