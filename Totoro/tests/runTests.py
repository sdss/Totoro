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

from __future__ import division, print_function

import sys
from subprocess import PIPE, Popen

import nose
from Totoro.core.colourPrint import _color_text


if __name__ == '__main__':
    """Runs Totoro test.

    IMPORTANT: nose should be called by using this script and not by just
    writing "nosetests". Otherwise, the production DB may be used, which may
    (or may not) break something.
    """

    if '-v' not in sys.argv:
        sys.argv.append('-v')
    if '--nologcapture' not in sys.argv:
        sys.argv.append('--nologcapture')
    if '--no-restore' not in sys.argv:
        # We try to run the script that restores the test DB. If it fails we
        # still run the test suite.
        sys.stdout.write('Restoring test DB ... ')
        sys.stdout.flush()
        command = Popen('restoreTestDB.sh', stdout=PIPE, stderr=PIPE, shell=True, close_fds=True)
        out, err = command.communicate()

        if err.decode('utf-8') != '':
            print(_color_text('FAILED', 'red'))
        else:
            print('ok')

    # We wait to import setDefaulProfile until now because otherwise
    # restoreTestDB would fail when trying to remove the DB, as there would be
    # an open connection.
    from Totoro.db import setDefaulProfile
    setDefaulProfile('test')

    nose.main()
