# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings
warnings.filterwarnings('ignore', module='astropy.time.core')

from .core.setupTools import generateVersion
from .exceptions import TotoroError

from .core.logger import initLog
log = initLog

__minimum_numpy_version__ = '1.5.0'
__minimum_astropy_version__ = '0.3.0'
__minimum_sqlalchemy_version__ = '0.9.0b1'

try:
    from .version import version as __version__
except:
    __version__ = ''

try:
    from .plateDB.dataModel import db, Base, session
except:
    raise TotoroError('Impossible to connect to plateDB.')

from .scheduler.scheduler import BaseScheduler, Planner

__ALL__ = ['db', 'Base', 'session', 'log', 'Scheduler']
