#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Based on the astropy installer.

from setuptools import setup
import os

import Totoro

NAME = 'Totoro'

# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
VERSION = '0.1.dev'

# Indicates if this version is a release version
RELEASE = 'dev' not in VERSION

# if not RELEASE:
#     VERSION += get_git_devstr(False)

DOWNLOAD_BASE_URL = 'https://github.com/albireox/Totoro'

# Freeze build information in version.py
Totoro.generateVersion(NAME, VERSION, RELEASE)

# Creates the package list
packages = []
for walkList in os.walk(NAME):
      packages.append(walkList[0])

setup_requires = ['numpy>=' + Totoro.__minimum_numpy_version__]
setup_requires += ['astropy>=' + Totoro.__minimum_astropy_version__]
setup_requires += ['sqlalchemy>=' + Totoro.__minimum_sqlalchemy_version__]
install_requires = setup_requires

setup(name=NAME,
      version=VERSION,
      packages=packages,
      description='A scheduling tool for MaNGA',
      requires=['numpy', 'astropy', 'sqlalchemy'],
      setup_requires=setup_requires,
      install_requires=install_requires,
      provides=[NAME],
      author='Jose R. Sanchez-Gallego',
      author_email='j.sanchezgallego@uky.edu',
      license='BSD',
      url='http://github.com/albireox/Totoro',
      zip_safe=False,
      use_2to3=True)
