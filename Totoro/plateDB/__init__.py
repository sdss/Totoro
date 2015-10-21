#!/usr/bin/env python
# encoding: utf-8
"""
__init__.py

Created by José Sánchez-Gallego on 30 Oct 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

import warnings
from sqlalchemy import exc as sa_exc

warnings.simplefilter('ignore', category=sa_exc.SAWarning)
warnings.filterwarnings('ignore', 'Predicate of partial index')

__all__ = ['Connect']

from .connection import DatabaseConnection
db = DatabaseConnection()
session = db.session
