#!/usr/bin/env python
# encoding: utf-8
"""
mangaPluggingRequest.py

Created by José Sánchez-Gallego on 15 Nov 2014.
Licensed under a 3-clause BSD license.

Revision history:
    15 Nov 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import os
import sys
import argparse
from astropy import time
from sdss.internal.manga.Totoro.exceptions import TotoroError
from sdss.internal.manga.Totoro.scheduler import Plugger
from sdss.internal.manga.Totoro.scheduler.observingPlan import ObservingPlan
from sdss.internal.manga.Totoro import log


def mangaPluggingRequest(mjd=None):
    """Creates a `Plugger` instance for the requested mjd and returns the
    output. If `mjd=None`, the current MJD will be used."""

    if mjd is None:
        tt = time.Time.now()
    else:
        tt = time.Time(int(mjd), format='mjd')

    printMJD = int(tt.jd) - 2400000

    log.info('creating plugging request for MJD={0:d}'.format(printMJD))
    observingPlan = ObservingPlan()
    jd0, jd1 = observingPlan.getJD(jd=int(tt.jd))

    if jd0 is None or jd1 is None:
        raise TotoroError('not MaNGA observations scheduled for '
                          'MJD={0:d}'.format(printMJD))

    plugger = Plugger(jd0, jd1)
    plugger.getOutput()


def main(argv=None):

    parser = argparse.ArgumentParser(description=__doc__,
                                     prog=os.path.basename(sys.argv[0]))
    parser.add_argument('-m', '--mjd', metavar='MJD', type=int, default=None,
                        help='The MJD for which the plugging is requested.')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose',
                        help='Print lots of extra output.')

    args = parser.parse_args()

    if args.verbose is True:
        log.sh.setLevel('DEBUG')

    mangaPluggingRequest(mjd=args.mjd)

    return


if __name__ == "__main__":
    main()
