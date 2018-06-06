#!/usr/bin/env python
# encoding: utf-8
#
# ancillaryPlatesPrio.py
#
# Created by José Sánchez-Gallego on 25 Apr 2016.
# Licensed under a 3-clause BSD license.
#
# Revision history:
#    25 Apr 2016 J. Sánchez-Gallego
#       Initial version

from __future__ import division, print_function

import numpy as np
from astropy import table

from Gohan.utils import getAllDrilledTargets, sdssMaskBits
from Totoro import getAll


ancillaryBits = sdssMaskBits[(sdssMaskBits['flag'] == 'MANGA_TARGET3') & (sdssMaskBits['bit'] > 0)]

allPlates = getAll()

targets = table.Table(getAllDrilledTargets())


def sortAncillaryPrograms(do_print=False):
    """Sorts the ancillary programs by number of observed targets."""

    ancillaries = targets[targets['manga_target3'] > 0]

    drilledPerProgram = []
    for mask in ancillaryBits:
        targetsInProgram = targets['manga_target3'] & 1 << mask['bit']
        drilledPerProgram.append(np.sum(targetsInProgram > 0))

    completePlates = [plate for plate in allPlates if plate.isComplete]
    plateids = [plate.plate_id for plate in completePlates]

    observedAncillaries = ancillaries[np.in1d(ancillaries['plateid'], plateids)]

    observedPerProgram = []
    for mask in ancillaryBits:
        targetsInProgram = (observedAncillaries['manga_target3'] & 1 << mask['bit'])
        observedPerProgram.append(np.sum(targetsInProgram > 0))

    ancillaryTable = table.Table(
        None, names=['label', 'bit', 'nDrilled', 'nObserved'], dtype=['S20', int, int, int])

    for ii, mask in enumerate(ancillaryBits):
        ancillaryTable.add_row((mask['label'], mask['bit'], drilledPerProgram[ii],
                                observedPerProgram[ii]))

    ancillaryTable.sort(['nObserved'])
    ancillaryTable.pprint()

    return ancillaryTable


def ancillaryPlatesPrio(sortedAncillary):
    """Returns a list of incomplete plates with ancillaries."""

    # Ancillary programs with fewer than 5 targets observed.
    ancillary5 = sortedAncillary[sortedAncillary['nObserved'] < 5]

    incompletePlates = [
        plate for plate in allPlates if not plate.isComplete and plate.priority > 1
    ]

    validPlates = []
    for plate in incompletePlates:
        plateTargets = targets[targets['plateid'] == plate.plate_id]
        for bit in ancillary5['bit']:
            if np.any(plateTargets['manga_target3'] & 1 << bit):
                validPlates.append(plate)
                break

    prioTable = table.Table(
        None,
        names=['plate_id', 'nTargets', 'labels', 'priority', 'location', 'status'],
        dtype=[int, int, 'S100', int, 'S20', 'S20'])

    for plate in validPlates:
        plateTargets = targets[targets['plateid'] == plate.plate_id]
        nTargets = 0
        labels = []
        for mask in ancillary5:
            bitTargets = plateTargets['manga_target3'] & 1 << mask['bit']
            nTargetsBit = np.sum(bitTargets > 0)
            nTargets += nTargetsBit
            if nTargetsBit > 0:
                labels.append('{0} ({1})'.format(mask['label'], nTargetsBit))

        prioTable.add_row((plate.plate_id, nTargets, ', '.join(labels), plate.priority,
                           plate.getLocation(), plate.statuses[0].label))

    print()

    prioTable.sort(['nTargets'])
    prioTable.reverse()
    prioTable.pprint(max_lines=100)


if __name__ == '__main__':
    sortedAncillary = sortAncillaryPrograms(do_print=True)
    ancillaryPlatesPrio(sortedAncillary)
