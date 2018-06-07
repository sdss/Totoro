#!/usr/bin/env python
# encoding: utf-8
"""

methods.py

Created by José Sánchez-Gallego on 19 Nov 2015.
Licensed under a 3-clause BSD license.

Revision history:
    19 Nov 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division, print_function

from sqlalchemy.orm.session import Session


def addFunctionsPlateDB(plateDB):

    def scienceExposuresPlugging(self):
        """Returns a list of science exposures for a plugging."""
        session = Session.object_session(self)
        exps = session.query(plateDB.Exposure).join(
            plateDB.Observation,
            plateDB.ExposureFlavor).filter(plateDB.Observation.plugging_pk == self.pk).filter(
                plateDB.ExposureFlavor.label == 'Science').all()
        return exps

    def scienceExposuresPlate(self):
        """Returns a list of science exposures for a plate."""
        exps = []
        for plugging in self.pluggings:
            exps += plugging.scienceExposures()
        return exps

    def mjd(self):
        """Returns the *SDSS* MJD.

        See line ~140 (the mjd4Gang function) here for notes on this value.
        https://svn.sdss.org/deprecated/operations/iop/trunk/etc/iopUtils.tcl
        """

        return int(float(self.start_time) / 86400.0 + 0.3)

    plateDB.Plugging.scienceExposures = scienceExposuresPlugging
    plateDB.Plate.scienceExposures = scienceExposuresPlate
    plateDB.Exposure.mjd = mjd
