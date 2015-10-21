#!/usr/bin/env python
# encoding: utf-8
"""
tabularOutput.py

Created by José Sánchez-Gallego on 5 Jul 2014.
Licensed under a 3-clause BSD license.

Revision history:
    5 Jul 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from astropy import table
import numpy as np
from sdss.internal.manga.Totoro.core.colourPrint import _color_text
import sys
import re

ansi_escape = re.compile(r'\x1b[^m]*m')


def printTabularOutput(plates, **kwargs):
    """Returns a series of tables with information about the schedule."""

    observedOutPuts = [ObservedPlateOutput(plate) for plate in plates]
    cartOrder = np.argsort([plate.getCartridgeNumber() for plate in plates])

    for idx in cartOrder:
        print()
        pText = observedOutPuts[idx].pprint()
        for line in pText:
            if sys.stdout.isatty():
                print(line)
            else:
                print(ansi_escape.sub('', line))

    return


def getTabularOutput(plates, **kwargs):
    """Returns a series of tables with information about the schedule."""

    observedOutPuts = [ObservedPlateOutput(plate) for plate in plates]
    cartOrder = np.argsort([plate.getCartridgeNumber() for plate in plates])

    text = ''

    for idx in cartOrder:
        text += '\n'
        pText = observedOutPuts[idx].pprint()
        for line in pText:
            text += line
            text += '\n'

    return text


class ObservedPlateOutput(dict):

    def __init__(self, plate, **kwargs):

        self.sets = []
        self.plate = plate

        super(ObservedPlateOutput, self).__init__(
            {'plate': ObservedTable(), 'sets': self.sets})

        self._createPlateColumns()
        self._createSets()

    def _createPlateColumns(self):

        self['plate'].addColumn(
            ObservedField([self.plate.plate_id], 'plate_id', int))

        isPlateComplete = self.plate.isComplete

        if isPlateComplete:
            text = 'YES'
        else:
            text = 'NO'

        self['plate'].addColumn(ObservedField([text], 'complete', 'S10'))

        self['plate'].addColumn(
            ObservedField([self.plate.getCartridgeNumber()], 'cart #', int))

        sn2Values = self.plate.getCumulatedSN2(includeIncomplete=True)
        self['plate'].addColumn(ObservedField([sn2Values[0]], 'b1', float))
        self['plate'].addColumn(ObservedField([sn2Values[1]], 'b2', float))
        self['plate'].addColumn(ObservedField([sn2Values[2]], 'r1', float))
        self['plate'].addColumn(ObservedField([sn2Values[3]], 'r2', float))

        self['plate'].addColumn(
            ObservedField(
                [len(self.plate.getValidSets(includeIncomplete=True))],
                '# sets', int))

        self['plate'].addColumn(
            ObservedField(
                [len(self.plate.getValidExposures())], '# exps', int))

        utRange = self.plate.getUTVisibilityWindow()
        self['plate'].addColumn(
            ObservedField([utRange[0]], 'minUT', 'S10'))
        self['plate'].addColumn(
            ObservedField([utRange[1]], 'maxUT', 'S10'))

    def _createSets(self):

        validSets = np.array(self.plate.getValidSets(includeIncomplete=True))
        incomplete = np.array(
            [set.getQuality() == 'Incomplete' for set in validSets])

        nn = 1
        for set in validSets[np.where(incomplete)]:
            sTable = ObservedSetOutput(set, nn)
            self.sets.append(sTable)
            nn += 1

        for set in validSets[np.where(incomplete == np.bool(False))]:
            sTable = ObservedSetOutput(set, nn)
            self.sets.append(sTable)
            nn += 1

    def pprint(self, printSets=True):

        platePFormat = self['plate'].pformat(max_lines=1e6, max_width=1e6)
        for ii in range(2):
            platePFormat[ii] = _color_text(platePFormat[ii], 'red')

        platePFormat.append('\n')
        if printSets:
            for set in self.sets:
                setPprint = set.pprint(printExposures=True)
                for line in setPprint:
                    platePFormat.append('      ' + line)
                platePFormat.append('\n\n')

        return platePFormat


class ObservedSetOutput(dict):

    def __init__(self, set, nSet, **kwargs):

        self.exposures = None
        self.set = set
        self.nSet = nSet

        super(ObservedSetOutput, self).__init__(
            {'set': ObservedTable(), 'exposures': self.exposures})

        self._createSetColumns()
        self._createExposures()

    def _createSetColumns(self):

        self['set'].addColumn(ObservedField([self.nSet], 'set #', 'S10'))

        if self.set.getQuality() != 'Bad':
            self['set'].addColumn(
                ObservedField([self.set.getQuality()], 'status', 'S10'))

        self['set'].addColumn(
            ObservedField(
                [len(self.set.getValidExposures())], '# exps', 'S10'))

        self['set'].addColumn(
            ObservedField(
                [', '.join(self.set.getDitherPositions())],
                'ditherPos', 'S10'))

        sn2Values = self.set.getSN2Array()
        self['set'].addColumn(ObservedField([sn2Values[0]], 'b1', float))
        self['set'].addColumn(ObservedField([sn2Values[1]], 'b2', float))
        self['set'].addColumn(ObservedField([sn2Values[2]], 'r1', float))
        self['set'].addColumn(ObservedField([sn2Values[3]], 'r2', float))

        self['set'].addColumn(
            ObservedField([self.set.getAverageSeeing()], 'avgSeeing', float))

        if self.set.getQuality() != 'Incomplete':
            self['set'].addColumn(ObservedField([-999.], 'minUT', float))
            self['set'].addColumn(ObservedField([-999.], 'maxUT', float))
        else:
            utRange = self.set.getUTVisibilityWindow()
            self['set'].addColumn(ObservedField([utRange[0]], 'minUT', 'S10'))
            self['set'].addColumn(ObservedField([utRange[1]], 'maxUT', 'S10'))

    def _createExposures(self):

        eTable = ObservedExposureOutput(self.set.getValidExposures())
        self.exposures = eTable
        self['exposures'] = eTable

    def pprint(self, printExposures=True):

        setPFormat = self['set'].pformat(max_lines=1e6, max_width=1e6)
        for ii in range(2):
            setPFormat[ii] = _color_text(setPFormat[ii], 'blue') + \
                _color_text('', 'default')

        if not printExposures:
            return setPFormat
        else:
            setPFormat.append('\n')
            expPprint = self.exposures.pprint()
            for line in expPprint:
                setPFormat.append('     ' + line)

            return setPFormat


class ObservedExposureOutput(object):

    def __init__(self, exp, **kwargs):

        self.exposures = exp

        self.expTable = ObservedTable()

        self._createExposureColumns()

    def _createExposureColumns(self):

        ditherPosCol = ObservedField(
            [exp.ditherPosition for exp in self.exposures], 'ditherPos', 'S10')
        self.expTable.addColumn(ditherPosCol)

        sn2Values = [exp.getSN2Array() for exp in self.exposures]
        self.expTable.addColumn(ObservedField(zip(*sn2Values)[0], 'b1', float))
        self.expTable.addColumn(ObservedField(zip(*sn2Values)[1], 'b2', float))
        self.expTable.addColumn(ObservedField(zip(*sn2Values)[2], 'r1', float))
        self.expTable.addColumn(ObservedField(zip(*sn2Values)[3], 'r2', float))

        self.expTable.addColumn(
            ObservedField([exp.manga_seeing for exp in self.exposures],
                          'seeing', float))

        utRange = [exp.getUTObserved(format='str') for exp in self.exposures]
        self.expTable.addColumn(
            ObservedField(zip(*utRange)[0], 'minUT', 'S10'))
        self.expTable.addColumn(
            ObservedField(zip(*utRange)[1], 'maxUT', 'S10'))

    def pprint(self, show_name=True):

        expPFormat = self.expTable.pformat(
            max_lines=1e6, max_width=1e6, show_name=show_name)

        for ii in range(2):
            expPFormat[ii] = _color_text(expPFormat[ii], 'yellow') + \
                _color_text('', 'default')

        return expPFormat


class ObservedTable(table.Table):

    def __init__(self):

        super(ObservedTable, self).__init__(None)

    def addColumn(self, column):

        self.add_column(column)


class ObservedField(table.Column):

    def __init__(self, data, name, dtype, **kwargs):

        if not isinstance(data, (list, tuple, np.ndarray)):
            super(ObservedField, self).__init__(
                data=[data], name=[name], dtype=[dtype])
        else:
            super(ObservedField, self).__init__(
                data=data, name=name, dtype=dtype)
