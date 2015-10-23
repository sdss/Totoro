#!/usr/bin/env python
# encoding: utf-8
"""

encryptProfile.py

Created by José Sánchez-Gallego on 10 Oct 2015.
Licensed under a 3-clause BSD license.

Revision history:
    10 Oct 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from simplecrypt import encrypt
import sys
import configparser
import os


def encryptProfile(file):
    """Encrypts a profile file."""

    if 'SDSSCONNECT_PASSWORD' not in os.environ:
        raise RuntimeError('$SDSSCONNECT_PASSWORD not defined')

    passwd = os.environ['SDSSCONNECT_PASSWORD']

    splitPath = os.path.splitext(file)
    outputPath = splitPath[0] + '_enc' + splitPath[1]

    if os.path.exists(outputPath):
        raise RuntimeError('output file {0} exists. Remove the file first if '
                           'you want to overwrite it.'.format(outputPath))

    # Checks that the input file can be read by configparser
    try:
        config = configparser.ConfigParser()
        config.read(file)
    except Exception as ee:
        raise RuntimeError('the input file does not have a valid format. '
                           'Error: {0}'.format(ee))

    plainText = open(file, 'r').read()
    with open(outputPath, 'wb') as output:
        cipherText = encrypt(passwd, plainText)
        output.write(cipherText)

    print('Profile file {0} saved to {1}'.format(file, outputPath))


if __name__ == '__main__':
    encryptProfile(sys.argv[1])
