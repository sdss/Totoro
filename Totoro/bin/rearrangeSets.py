#!/usr/bin/env python
# encoding: utf-8
"""
rearrageSets.py

Created by José Sánchez-Gallego on 10 Dec 2014.
Licensed under a 3-clause BSD license.

Revision history:
    10 Dec 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division, print_function

import argparse
import os
import sys
import warnings


TotoroPath = os.path.realpath(os.path.join(os.path.dirname(__file__), '../../'))
if TotoroPath not in sys.path:
    sys.path.append(TotoroPath)


def rearrageSets(plateid, force=False, LST=None, **kwargs):
    """Triggers a set rearrangement for a certain plate."""

    from Totoro.dbclasses import plate

    pp = plate.fromPlateID(plateid)
    pp.rearrangeSets(force=force, LST=LST)


if __name__ == '__main__':

    warnings.simplefilter('always', DeprecationWarning)
    warnings.warn('this script is now deprecated an may be removed in a future '
                  'version of Totoro. Please, use "totoro rearrange" instead.',
                  DeprecationWarning)
    warnings.simplefilter('ignore', DeprecationWarning)

    parser = argparse.ArgumentParser(description=__doc__, prog=os.path.basename(sys.argv[0]))
    parser.add_argument(
        'plate', metavar='PLATE', type=int, help='the plates which sets will be rearranged')
    parser.add_argument(
        '-f',
        '--force',
        action='store_true',
        dest='force',
        help='forces the set rearrangement even if the number '
        'of exposures exceeds the limit.')
    parser.add_argument(
        '-o',
        '--outputflush',
        action='store_true',
        dest='outputflush',
        help='force sdtout writes to be immediately flushed '
        'used for unix redirection and file monitoring')
    parser.add_argument(
        '-n',
        '--nolst',
        action='store_false',
        dest='LST',
        help='does not use the current LST to rearrange the '
        'sets. If set, the incomplete sets will be moved '
        'at the beginning of the observing window, if '
        'possible.')
    parser.add_argument(
        '-l',
        '--lst',
        type=float,
        dest='LST',
        help='tries to rearrange sets so that at least an '
        'incomplete set is available after that LST.')
    parser.add_argument(
        '-v', '--verbose', action='store_true', dest='verbose', help='Print lots of extra output.')

    args = parser.parse_args()

    if args.outputflush:
        sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    rearrageSets(args.plate, force=args.force, LST=args.LST)
