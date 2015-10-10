#!/usr/bin/env python
# encoding: utf-8
"""

plateDB.py

Created by José Sánchez-Gallego on 10 Oct 2015.
Licensed under a 3-clause BSD license.

Revision history:
    10 Oct 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import os
import pickle


def construct_plateDB(Base):

    metadata_pickle_filename = "ModelClasses_apo_platedb.pickle"
    use_cache = 'MODEL_CLASSES_CACHE' in os.environ \
        and os.environ['MODEL_CLASSES_CACHE'] == 'YES'

    if use_cache:
        cache_path = os.path.join(os.path.expanduser('~'), '.sqlalchemy_cache')
        cached_metadata = None
        if os.path.exists(cache_path):
            try:
                cached_metadata = pickle.load(
                    open(os.path.join(cache_path, metadata_pickle_filename)))
            except IOError:
                # cache file not found
                pass
    else:
        cached_metadata = False

    class PlateDB(object):
        class ActivePlugging(Base):
            if cached_metadata:
                __table__ = cached_metadata.tables['platedb.active_plugging']
            else:
                __tablename__ = 'active_plugging'
                __table_args__ = {'autoload': True, 'schema': 'platedb'}

            def __repr__(self):
                return '<Active Plugging: plate={0} cartridge={1}>'.format(
                        self.plugging.plate.plate_id,
                        self.plugging.cartridge.number)

    return PlateDB
