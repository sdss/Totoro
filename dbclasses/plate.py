#!/usr/bin/env python
# encoding: utf-8
"""
plate.py

Created by José Sánchez-Gallego on 18 Apr 2014.
Licensed under a 3-clause BSD license.

Revision history:
    18 Apr 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from ..apoDB import TotoroDBConnection
from ..utils import mlhalimit, isPlateComplete
from ..exceptions import TotoroError, MultiplePlatePointings, PlateNotPlugged
from ..exceptions import NoMangaExposure
import warnings
from .. import log, config, dustMap, site
from .set import getPlateSets, Set
from .exposure import Exposure
from ..logic import setArrangement, updateSets, getValidSet
import numpy as np
from ..utils import getIntervalIntersectionLength, getAPOcomplete
from copy import deepcopy
import os


totoroDB = TotoroDBConnection()
plateDB = totoroDB.plateDB
mangaDB = totoroDB.mangaDB
session = totoroDB.Session()


class Plates(list):

    def __init__(self, inp, format='pk', **kwargs):

        if all([isinstance(ii, Plate) for ii in inp]):
            list.__init__(self, inp)
        else:
            list.__init__(self, [Plate(ii, format=format, **kwargs)
                                 for ii in inp])

    @staticmethod
    def getPlugged(onlyIncomplete=False, **kwargs):

        with session.begin(subtransactions=True):
            activePluggings = session.query(
                plateDB.ActivePlugging).join(
                    plateDB.Plugging,
                    plateDB.Plate,
                    plateDB.PlateToSurvey,
                    plateDB.Survey, plateDB.SurveyMode).filter(
                        plateDB.Survey.label == 'MaNGA',
                        plateDB.SurveyMode.label == 'MaNGA dither').order_by(
                            plateDB.Plate.plate_id).all()

        plates = [actPlug.plugging.plate_pk for actPlug in activePluggings]

        if onlyIncomplete:
            return Plates._getIncomplete(plates, **kwargs)
        else:
            return Plates(plates, **kwargs)

    @staticmethod
    def getAtAPO(onlyIncomplete=False, **kwargs):

        with session.begin(subtransactions=True):
            plates = session.query(plateDB.Plate).join(
                plateDB.PlateLocation).filter(
                    plateDB.PlateLocation.label == 'APO').join(
                        plateDB.PlateToSurvey).join(
                            plateDB.Survey, plateDB.SurveyMode).filter(
                                plateDB.Survey.label == 'MaNGA',
                                plateDB.SurveyMode.label ==
                                'MaNGA dither').order_by(
                                    plateDB.Plate.plate_id).all()

        plates = [plate.pk for plate in plates]

        if onlyIncomplete:
            return Plates._getIncomplete(plates, **kwargs)
        else:
            return Plates(plates, **kwargs)

    @staticmethod
    def getAll(onlyIncomplete=False, **kwargs):

        with session.begin(subtransactions=True):
            plates = session.query(plateDB.Plate).join(
                plateDB.PlateToSurvey, plateDB.Survey, plateDB.SurveyMode
                ).filter(plateDB.Survey.label == 'MaNGA',
                         plateDB.SurveyMode.label == 'MaNGA dither').order_by(
                    plateDB.Plate.plate_id).all()

        plates = [plate.pk for plate in plates]

        if onlyIncomplete:
            return Plates._getIncomplete(plates, **kwargs)
        else:
            return Plates(plates, **kwargs)

    @staticmethod
    def _getIncomplete(plates, **kwargs):

        totoroPlates = [Plate(plate, format='pk') for plate in plates]

        incompletePlates = []
        for totoroPlate in totoroPlates:
            if not isPlateComplete(totoroPlate):
                incompletePlates.append(totoroPlate)
        return incompletePlates


class Plate(plateDB.Plate):

    def __new__(cls, input=None, format='pk', **kwargs):

        if input is None:
            return plateDB.Plate.__new__(cls)

        base = cls.__bases__[0]

        with session.begin(subtransactions=True):
            instance = session.query(base).filter(
                eval('{0}.{1} == {2}'.format(base.__name__, format, input))
                ).one()

        instance.__class__ = cls

        return instance

    def __init__(self, input, format='pk', mock=False, silent=False,
                 updateSets=True, **kwargs):

        self._complete = None
        self.isMock = mock
        self.kwargs = kwargs

        if 'dust' in kwargs:
            self.dust = kwargs['dust']
        else:
            self.dust = None if dustMap is None else dustMap(self.ra, self.dec)

        self.mlhalimit = mlhalimit(self.dec)

        if not self.isMock:
            self.checkPlate()
            self.sets = getPlateSets(self.pk, format='pk', silent=silent,
                                     **kwargs)
            if updateSets:
                self.updateSets()

            if silent is False:
                log.debug('loaded plate with pk={0}, plateid={1}'.format(
                          self.pk, self.plate_id))

        else:
            self.sets = []

    def __repr__(self):
        return ('<Totoro Plate (plate_id={0}, pk={1}, completion={2:.2f})'
                .format(self.plate_id, self.pk, self.getPlateCompletion()))

    @classmethod
    def fromPlateID(cls, plateid, **kwargs):
        return Plate(plateid, format='plate_id', **kwargs)

    @classmethod
    def fromSets(cls, sets, silent=False, **kwargs):

        newPlate = plateDB.Plate.__new__(cls)
        newPlate.__init__(None, mock=True, **kwargs)
        newPlate.sets = sets

        if not silent:
            log.debug('created mock plate from sets pk={0}'.format(
                ', '.join(map(str, sets))))

        return newPlate

    @classmethod
    def createMockPlate(cls, ra=None, dec=None, silent=False, **kwargs):

        if ra is None or dec is None:
            raise TotoroError('ra and dec must be specified')

        mockPlate = plateDB.Plate.__new__(cls)
        mockPlate.__init__(None, mock=True, ra=ra, dec=dec, **kwargs)

        mockPlate.isMock = True

        if not silent:
            log.debug(
                'created mock plate with ra={0:.3f} and dec={0:.3f}'.format(
                    mockPlate.coords[0], mockPlate.coords[1]))

        return mockPlate

    def checkPlate(self):

        if not hasattr(self, 'pk') or not hasattr(self, 'plate_id'):
            raise AttributeError('Plate instance has no pk or plate_id.')

        if not self.isMaNGA:
            raise TotoroError('this is not a MaNGA plate!')

        nMaNGAExposures = len(self.getMangadbExposures())
        nScienceExposures = len(self.getScienceExposures())
        if nMaNGAExposures != nScienceExposures:
            warnings.warn('{0} plateDB.Exposures found for plate_id={1} '
                          'but only {2} mangaDB.Exposures'.format(
                              nScienceExposures, self.plate_id,
                              nMaNGAExposures), NoMangaExposure)

    @property
    def isMaNGA(self):

        for survey in self.surveys:
            if survey.label == 'MaNGA':
                return True

        return False

    def updateSets(self):
        result = updateSets(self)
        if result:
            self.update()

    def rearrangeSets(self, startDate=None, **kwargs):
        result = setArrangement.rearrangeSets(
            self, startDate=startDate, **kwargs)
        if result is True:
            self.update()

    @property
    def ra(self):
        return self.getCoordinates()[0]

    @property
    def dec(self):
        return self.getCoordinates()[1]

    @property
    def coords(self):
        """Alias for getCoordinates(). Provides backwards compatibility."""
        return self.getCoordinates()

    def getCoordinates(self):

        if 'ra' in self.kwargs and 'dec' in self.kwargs:
            if (self.kwargs['ra'] is not None and
                    self.kwargs['dec'] is not None):
                return np.array(
                    [self.kwargs['ra'], self.kwargs['dec']], np.float)
        else:

            if len(self.plate_pointings) > 1:
                warnings.warn('multiple plate pointings for plate_id={0:d}. '
                              'Using the first one.'.format(self.plate_id),
                              MultiplePlatePointings)

            return np.array(
                [self.plate_pointings[0].pointing.center_ra,
                 self.plate_pointings[0].pointing.center_dec], np.float)

    @property
    def isComplete(self):
        if self._complete is not None:
            return self._complete
        else:
            return isPlateComplete(self)

    def copy(self):
        return deepcopy(self)

    def getPlateCompletion(self, includeIncompleteSets=False):

        totalSN = self.getCumulatedSN2(
            includeIncomplete=includeIncompleteSets)

        completionRates = totalSN.copy()
        completionRates[0:2] /= config['SN2thresholds']['plateBlue']
        completionRates[2:] /= config['SN2thresholds']['plateRed']

        return np.min(completionRates)

    def getCumulatedSN2(self, includeIncomplete=False):

        validStatuses = ['Good', 'Excellent']
        if includeIncomplete:
            validStatuses.append('Incomplete')

        validSets = []
        for set in self.sets:
            if set.getQuality()[0] in validStatuses:
                validSets.append(set)

        if len(validSets) == 0:
            return np.array([0.0, 0.0, 0.0, 0.0])
        else:
            return np.sum([set.getSN2Array() for set in validSets], axis=0)

    def getActiveCartNumber(self):

        for plugging in self.pluggings:
            if len(plugging.activePlugging) > 0:
                return int(plugging.cartridge.number)

        raise PlateNotPlugged(
            'plate_id={0} is not currently plugged'.format(self.plate_id))

    def getMangadbExposures(self):
        """Returns a list of mangaDB.Exposure objects with all the exposures
        for this plate."""

        with session.begin(subtransactions=True):
            mangaExposures = session.query(totoroDB.mangaDB.Exposure).join(
                plateDB.Exposure,
                plateDB.Observation,
                plateDB.Plugging,
                plateDB.Plate).filter(
                    plateDB.Plate.pk == self.pk).all()

        return mangaExposures

    def getScienceExposures(self):
        """Returns a list of all plateDB.Exposure science exposures for
        the plate."""

        scienceExps = []

        for plugging in self.pluggings:
            scienceExps += plugging.scienceExposures()

        return scienceExps

    def getTotoroExposures(self):
        """Returns a list of Totoro dbclasses.Exposure instances for this
        Plate insance."""

        totoroExposures = []
        for set in self.sets:
            totoroExposures += set.totoroExposures

        return totoroExposures

    def getValidSets(self, includeIncomplete=False):

        validSets = []
        for set in self.sets:
            quality = set.getQuality()
            if quality in ['Good', 'Excellent', 'Poor']:
                validSets.append(set)
            elif includeIncomplete and quality == 'Incomplete':
                validSets.append(set)

        return validSets

    def update(self, **kwargs):
        """Reloads the plate."""

        kwargs.update(self.kwargs)

        newSelf = Plate(self.pk, format='pk', **kwargs)
        self = newSelf

        log.debug('plate with plate_id={0} has been reloaded'.format(
                  self.plate_id))

    def getValidExposures(self):
        """Returns all valid exposures, even if they belong to an incomplete
        or bad set."""

        validExposures = []

        for set in self.sets:
            for exp in set.totoroExposures:
                if exp.valid is True:
                    validExposures.append(exp)

        return validExposures

    def getLSTRange(self):

        ha0 = -self.mlhalimit
        ha1 = self.mlhalimit

        lst0 = (ha0 + self.coords[0]) % 360. / 15
        lst1 = (ha1 + self.coords[0]) % 360. / 15

        return (lst0, lst1)

    def getUTVisibilityWindow(self, format='str'):

        lst0, lst1 = self.getLSTRange()

        ut0 = site.localTime(lst0, utc=True, returntype='datetime')
        ut1 = site.localTime(lst1, utc=True, returntype='datetime')

        if format == 'str':
            return ('{0:%H:%M}'.format(ut0), '{0:%H:%M}'.format(ut1))
        else:
            return (ut0, ut1)

    def isVisible(self, LST0, LST1, minLength=None):
        """Returns True if the plate is visible in a range of LST."""

        if minLength is None:
            minLength = config['exposure']['exposureTime']

        lstIntersectionLength = getIntervalIntersectionLength(
            self.getLSTRange(), (LST0, LST1), wrapAt=24.)

        secIntersectionLength = lstIntersectionLength * 3600.

        return True if secIntersectionLength >= minLength else False

    def addMockExposure(self, exposure=None, startTime=None, set=None,
                        expTime=None, silent=False, **kwargs):
        """Creates a mock expusure in the best possible way."""

        ra, dec = self.coords

        if exposure is None:
            exposure = Exposure.createMockExposure(
                startTime=startTime, expTime=expTime, ra=ra, dec=dec,
                silent=silent, **kwargs)

        if exposure.isValid(silent=silent)[0] is False:
            if not silent:
                log.debug('mock exposure is invalid.')
            return False

        validSet = getValidSet(exposure, self)
        added = False
        for set in self.sets:
            if set == validSet:
                set.totoroExposures.append(exposure)
                added = True
        if not added:
            self.sets.append(Set.fromExposures([exposure], sient=silent))

        return exposure

    def getIncompleteSet(self, startTime, expTime):
        """Returns incomplete sets that are valid for an exposure starting at
        JD=startTime."""

        LST0 = site.localSiderialTime(startTime)
        LST1 = (LST0 + expTime / 60.) % 24.

        incompleteSets = [set for set in self.sets if not set.complete]

        if len(incompleteSets) == 0:
            return None

        for set in incompleteSets:

            lstRange = set.getLSTRange()
            intersectionLength = getIntervalIntersectionLength(
                (LST0, LST1), lstRange)

            if intersectionLength * 3600 > expTime:
                return set

        return None

    def getLastExposure(self):
        """Returns the last exposure taken."""

        exposures = []
        for set in self.sets:
            for exp in set.totoroExposures:
                if exp.valid:
                    exposures.append(exp)

        startTime = [exp.start_time for exp in exposures]
        order = np.argsort(startTime)

        return exposures[order[-1]]

    def getPriority(self):
        return self.plate_pointings[0].priority

    @property
    def isPlugged(self):
        try:
            self.getActiveCartNumber()
            return True
        except:
            return False

    def getAPOcomplete(self):
        return getAPOcomplete(self)

    def getMangaTileID(self):
        if len(self.design.inputs) == 0:
            return None
        else:
            for input in self.design.inputs:
                basename = os.path.basename(input.filepath)
                if 'mangaScience' not in basename:
                    continue
                noExt = os.path.splitext(basename)[0]
                return int(noExt.split('_')[1])
            return None
