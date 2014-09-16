#!/usr/bin/env python
# encoding: utf-8
"""
testPlugger.py

Created by José Sánchez-Gallego on 12 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    12 Aug 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from Totoro.scheduler import Plugger


def testPlugger():

    pl = Plugger(date=2456893.0)
    pl.getOutput()


if __name__ == '__main__':
    testPlugger()
