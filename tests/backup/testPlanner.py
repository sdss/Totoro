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
from Totoro.scheduler import Planner


def testPlanner():

    drilledLocIDs = [4863, 4867, 4869, 4870, 4871, 4872, 5668, 5669, 5670,
                     5671, 5672, 5673, 5676, 5737, 5738, 5740, 5741, 5742,
                     5743, 5744, 5745, 5746, 5747, 5748, 5749, 5751, 5771,
                     5772, 5773, 5774, 5783, 5784, 5805, 5806, 5807, 5808,
                     5809, 5802, 5845]

    pp = Planner(startDate=2456918.0, endDate=2456955.004861)

    for locID in drilledLocIDs:
        pp.fields.removeField(locID)

    pp.scheduleTimelines()

    return pp


if __name__ == '__main__':
    testPlanner()
