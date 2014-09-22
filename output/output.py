#!/usr/bin/env python
# encoding: utf-8
"""
output.py

Created by José Sánchez-Gallego on 22 Jul 2014.
Licensed under a 3-clause BSD license.

Revision history:
    22 Jul 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from collections import OrderedDict
from sdss.internal.manga.Totoro import scheduler
import numpy as np
import yaml
from astropy import time


__ALL__ = ['getNightlyOutput']


def formatValue(value):

    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def getNightlyOutput(input=None, format='dict', **kwargs):

    if input is None:
        nightly = scheduler.Nightly(**kwargs)
    else:
        nightly = input

    # Determines the MJD for which the windows are calculated
    if len(nightly.observingBlocks) == 0:
        mjd = None
    else:
        jd = np.mean([nightly.observingBlocks['JD0'],
                      nightly.observingBlocks['JD1']])
        mjd = time.Time(jd, format='jd', scale='tai').mjd

    output = OrderedDict()

    output['mangaStart'] = float(nightly.startDate) \
        if nightly.startDate is not None else None
    output['mangaEnd'] = float(nightly.endDate) \
        if nightly.endDate is not None else None
    output['currentTime'] = float(nightly.currentDate)

    # For future implementation
    output['schedule'] = None

    output['plates'] = OrderedDict()

    for plate in nightly.plates:

        output['plates'][plate.plate_id] = OrderedDict()
        thisPlate = output['plates'][plate.plate_id]

        thisPlate['cartridge'] = plate.getActiveCartNumber()
        thisPlate['complete'] = plate.isComplete
        thisPlate['completionPercentage'] = plate.getPlateCompletion() * 100.
        thisPlate['HARange'] = formatValue(plate.getHARange(intersect=False))
        thisPlate['UTVisibilityWindow'] = formatValue(
            plate.getUTRange(mjd=mjd, intersect=False, returnType='datetime'))
        thisPlate['SN2'] = formatValue(plate.getCumulatedSN2())

        thisPlate['sets'] = OrderedDict()

        for set in plate.sets:

            thisPlate['sets'][set.pk] = OrderedDict()
            thisSet = thisPlate['sets'][set.pk]

            thisSet['complete'] = set.complete

            thisSet['averageSeeing'] = formatValue(set.getAverageSeeing())
            thisSet['SN2'] = formatValue(set.getSN2Array())

            thisSet['SN2Range'] = formatValue(set.getSN2Range())
            thisSet['seeingRange'] = formatValue(set.getSeeingRange())
            thisSet['UTRange'] = set.getUTVisibilityWindow(
                mjd=mjd, intersect=False, returnType='datetime')
            thisSet['HARange'] = formatValue(set.getHARange(intersect=False))
            thisSet['missingDithers'] = formatValue(
                set.getMissingDitherPositions())

            thisSet['exposures'] = OrderedDict()

            for exposure in set.totoroExposures:

                thisSet['exposures'][exposure.pk] = OrderedDict()
                thisExposure = thisSet['exposures'][exposure.pk]

                thisExposure['valid'] = exposure.valid
                thisExposure['ditherPosition'] = exposure.ditherPosition
                thisExposure['obsHARange'] = formatValue(exposure.getHA())
                thisExposure['obsJDRange'] = formatValue(exposure.getJD())

                thisExposure['seeing'] = exposure.seeing
                thisExposure['SN2'] = formatValue(exposure.getSN2Array())

        if thisPlate['complete']:
            thisPlate['incompleteSets'] = []
        else:
            thisPlate['incompleteSets'] = OrderedDict(
                [(set, thisPlate['sets'][set])
                 for set in thisPlate['sets']
                 if thisPlate['sets'][set]['complete'] is False])

    if format == 'yaml':
        return yaml.dump(output, default_flow_style=False)

    return output
