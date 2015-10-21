#!/usr/bin/env python
# encoding: utf-8
"""
testAPOcomplete.py

Created by José Sánchez-Gallego on 14 Jul 2014.
Licensed under a 3-clause BSD license.

Revision history:
    14 Jul 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function

from Totoro.utils import getAPOcomplete, isPlateComplete

print(isPlateComplete(7815))
print(getAPOcomplete(7815))
