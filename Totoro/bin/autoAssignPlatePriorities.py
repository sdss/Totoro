#!/usr/bin/env python
# encoding: utf-8
"""
autoAssignPlatePriorities.py

Created by José Sánchez-Gallego on 11 Mar 2015.
Licensed under a 3-clause BSD license.

Revision history:
    11 Mar 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from Totoro.dbclasses import getAll
from Totoro import TotoroDBConnection
import numpy as np
import argparse
import os
import sys


db = TotoroDBConnection()

points = np.array([[-15, 4],
                   [-10, 4],
                   [-5, 4],
                   [0, 5],
                   [5, 5],
                   [10, 5],
                   [15, 5],
                   [20, 6],
                   [25, 4],
                   [30, 3],
                   [32, 3],
                   [35, 3],
                   [40, 4],
                   [43, 7],
                   [45, 7],
                   [50, 6],
                   [55, 5],
                   [60, 5],
                   [65, 5],
                   [70, 5]])

pp = np.polyfit(points[:, 0], points[:, 1], 10)


def calculatePriority(dec):
    """Calculates the priority for a given declination."""

    # Doesn't do anything for plates that are likely in the SGC.
    if dec < 22:
        return 5

    priority = int(np.round(np.polyval(pp, dec)))

    # Hard limits
    if priority > 7:
        priority = 7
    if priority < 3:
        priority = 3

    return priority


def autoAssignPlatePriorities(dryRun=False):
    """Assigns priorities to plates in the DB based on their declination."""

    allPlates = getAll(onlyIncomplete=True, rejectSpecial=True)

    for plate in allPlates:
        if plate.priority != 5:
            # Skips plates with manually assigned priority.
            continue
        newPriority = calculatePriority(plate.dec)
        if not dryRun:
            with db.session.begin():
                plate.plate_pointings[0].priority = newPriority
        print('plate_id={0} (Dec={1:.2f}) assigned priority={2}'.format(
              plate.plate_id, plate.dec, newPriority))

    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                                     prog=os.path.basename(sys.argv[0]))
    parser.add_argument('--dry-run', action='store_true', dest='dryrun',
                        help='Only tests changes.')

    args = parser.parse_args()

    autoAssignPlatePriorities(dryRun=args.dryrun)
