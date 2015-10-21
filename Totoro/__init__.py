# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings
warnings.filterwarnings('ignore', module='astropy.time.core')
warnings.filterwarnings(
    'ignore', 'Module argparse was already imported')

from .core.logger import initLog
log = initLog()

__minimum_numpy_version__ = '1.5.0'
__minimum_astropy_version__ = '0.3.0'
__minimum_sqlalchemy_version__ = '0.9.0b1'

try:
    from .version import version as __version__
except:
    __version__ = ''

from .scheduler import Planner, Field, Fields, Set, Exposure
