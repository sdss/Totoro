# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings
from .exceptions import DustMapWarning, TotoroError
# warnings.filterwarnings('ignore', module='astropy.time.core')
warnings.filterwarnings('ignore', 'Module argparse was already imported')

from .readPath import readPath
__DEFAULT_CONFIG_FILE__ = readPath('+defaults.yaml')
__TOTORO_CONFIG_PATH__ = readPath('~/.totoro/totoro.yaml')

# Reads the configuration file
from .core.configuration import getConfiguration
config = getConfiguration()

# Creates the custom logging system
from .core.logger import initLog
log = initLog()
log.debug('Logging starts now.')
log.debug('Configuration file has been loaded.')

try:
    from sdss.manga import DustMap
    dustMap = DustMap()
except (ImportError, ValueError):
    warnings.warn('no dust map found. No Galactic extinction '
                  'will be applied', DustMapWarning)
    dustMap = None
except:
    raise TotoroError('something went wrong while importing the dust map.')

from .utils import createSite
site = createSite()

from .apoDB import TotoroDBConnection
