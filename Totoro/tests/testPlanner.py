#!/usr/bin/env python
# encoding: utf-8
"""
testPlanner.py

Created by José Sánchez-Gallego on 11 Dec 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

from __future__ import division
from __future__ import print_function
from Totoro import Planner

from Totoro.core.logger import initLog


def testPlanner():

    initLog()
    pp = Planner()
    pp.simulate()

    return pp

if __name__ == '__main__':
    testPlanner()
