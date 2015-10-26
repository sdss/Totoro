#!/usr/bin/env python
# encoding: utf-8
"""
runTests.py

Created by José Sánchez-Gallego on 7 May 2015.
Licensed under a 3-clause BSD license.

Revision history:
    7 May 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from Totoro.db import setDefaulProfile
import nose
import sys


if __name__ == '__main__':
    """Runs Totoro test.

    IMPORTANT: nose should be called by using this script and not by just
    writing "nosetests". Otherwise, the production DB may be used, which may
    (or may not) break something.
    """

    setDefaulProfile('test')

    if '-v' not in sys.argv:
        sys.argv.append('-v')
    if '--nologcapture' not in sys.argv:
        sys.argv.append('--nologcapture')

    nose.main()
