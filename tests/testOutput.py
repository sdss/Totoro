#!/usr/bin/env python
# encoding: utf-8
"""
testOutput.py

Created by José Sánchez-Gallego on 23 Jul 2014.
Licensed under a 3-clause BSD license.

Revision history:
    23 Jul 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from Totoro.scheduler import Nightly
from Totoro.dbclasses import Plate
import cPickle


def testOutput():

    startDate = 2456843.1
    plates = [Plate(7815, format='plate_id',
                    rearrageExposures=True),
              Plate(7808, format='plate_id',
                    rearrageExposures=True)]

    nightly = Nightly(startDate=startDate, plates=plates)

    output = nightly.printTabularOutput()

    # blob = open('testOutput.pckl', 'w')
    # cPickle.dump(output, blob)
    # blob.close()

    return


if __name__ == '__main__':
    testOutput()
