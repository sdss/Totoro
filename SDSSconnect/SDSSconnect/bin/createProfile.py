#!/usr/bin/env python
# encoding: utf-8
"""

createProfile.py

Created by José Sánchez-Gallego on 10 Oct 2015.
Licensed under a 3-clause BSD license.

Revision history:
    10 Oct 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from simplecrypt import encrypt, decrypt
import StringIO
import getpass
import configparser
import os


def createProfileInteractive():
    """Interactively creates a profile."""

    profilesPath = os.path.join(os.path.expanduser('~'),
                                '.sdssconnect', 'profiles.ini')

    if 'SDSSCONNECT_PASSWORD' not in os.environ:
        raise RuntimeError('$SDSSCONNECT_PASSWORD not defined')

    passwd = os.environ['SDSSCONNECT_PASSWORD']

    profile = raw_input('Profile name: ').lower()

    config = configparser.ConfigParser()

    # If the profile file exists, reads it and checks if the profile we want to
    # add already exists.
    if os.path.exists(profilesPath):
        config.readfp(StringIO.StringIO(decrypt(passwd, open(profilesPath, 'r')
                      .read()).decode('utf8')))
        if config.has_section(profile):
            print('WARNING: profile exists. Overwriting it.')
        else:
            config.add_section(profile)
    else:
        config.add_section(profile)

    # Reads information
    database = raw_input('Database name: ')
    host = raw_input('Host: ')
    port = raw_input('Port: ')
    username = raw_input('Username: ')
    dbPasswd = getpass.getpass('Password: ')

    config.set(profile, 'user', username)
    config.set(profile, 'password', dbPasswd)
    config.set(profile, 'host', host)
    config.set(profile, 'port', port)
    config.set(profile, 'database', database)

    with open(profilesPath, 'wb') as output:
        buf = StringIO.StringIO()
        config.write(buf)
        ciphertext = encrypt(passwd, buf.getvalue())
        output.write(ciphertext)

    print('Profile {0} saved to {1}'.format(profile, profilesPath))


if __name__ == '__main__':
    createProfileInteractive()
