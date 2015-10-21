#!/usr/bin/env python
# encoding: utf-8
"""
testPlanner.py

Created by José Sánchez-Gallego on 11 Dec 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

from Totoro import Planner
from astropy import time

# Creates an observing plan in UTC
# ss = Scheduler(mode='planner')
# ss.observingPlan.save('nightly.D_UTC.txt')

starTime = time.Time(2456898.68, format='jd', scale='utc')
endTime = time.Time(2456916.7, format='jd', scale='utc')
ss = Planner(startTime=starTime, endTime=endTime)
# print ss.observingPlan
