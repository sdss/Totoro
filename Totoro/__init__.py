#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2018-06-07
# @Filename: __init__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# @Copyright: José Sánchez-Gallego

from __future__ import absolute_import

import warnings

from .exceptions import DustMapWarning, TotoroError
from .readPath import readPath


# warnings.filterwarnings('ignore', module='astropy.time.core')
warnings.filterwarnings('ignore', 'Module argparse was already imported')
warnings.filterwarnings(
    'ignore', 'Skipped unsupported reflection of ' + 'expression-based index q3c_field_idx')

__DEFAULT_CONFIG_FILE__ = readPath('+defaults.yaml')
__TOTORO_CONFIG_PATH__ = readPath('~/.totoro/totoro.yaml')

from .core.configuration import getConfiguration  # noqa
config = getConfiguration()

# Creates the custom logging system
from .core.logger import log  # noqa
log.start_file_logger(config['logging']['logFilePath'])
log.debug('Logging starts now.')
log.debug('Configuration file has been loaded.')

try:
    from .utils.dust_map import DustMap
    dustMap = DustMap()
except (ImportError, ValueError):
    log.warning('no dust map found. No Galactic extinction will be applied', DustMapWarning)
    dustMap = None
except Exception:
    raise TotoroError('something went wrong while importing the dust map.')

from .utils.site import Site  # noqa
site = Site()

__version__ = '2.1.2'
