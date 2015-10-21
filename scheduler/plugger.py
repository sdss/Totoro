#!/usr/bin/env python
# encoding: utf-8
"""
pluger.py

Created by José Sánchez-Gallego on 20 Oct 2014.
Licensed under a 3-clause BSD license.

Revision history:
    20 Oct 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.internal.manga.Totoro import log, config, TotoroDBConnection
from sdss.internal.manga.Totoro.scheduler.timeline import Timeline
from sdss.internal.manga.Totoro import exceptions
from sdss.internal.manga.Totoro import utils
from collections import OrderedDict
import warnings
import numpy as np


db = TotoroDBConnection()


def getAvailableCart(carts):
    """Gets a cart without a plugged plate. If none is available, returns
    None. If a cart is available, returns a tuple (cartNumber, status) with
    status a string indicating why the cart is available. The input argument,
    carts, must be a list of cart numbers."""

    from sdss.internal.manga.Totoro.dbclasses import Plate

    session = db.Session()
    with session.begin():
        activePluggings = session.query(db.plateDB.ActivePlugging).all()

    usedCarts = [aP.plugging.cartridge.number for aP in activePluggings
                 if aP.plugging.cartridge.number in carts]
    unusedCarts = [cart for cart in carts if cart not in usedCarts]

    if len(unusedCarts) > 0:
        return (unusedCarts[0], 'empty')

    for aP in activePluggings:

        cartNumber = aP.plugging.cartridge.number

        if cartNumber not in carts:
            continue

        isMaNGAPlate = (aP.plugging.plate.currentSurveyMode is not None and
                        'MaNGA' in aP.plugging.plate.currentSurveyMode.label)

        if not isMaNGAPlate:
            return (cartNumber, 'noMaNGAPlate')

        totoroPlate = Plate(aP.plugging.plate.plate_id, format='plate_id')

        if utils.isPlateComplete(totoroPlate):
            return (cartNumber, 'complete')

        elif totoroPlate.getPlateCompletion() == 0:
            return (cartNumber, 'notStarted')

    return None


def getCartPriority(carts):
    """Returns the priority of a list of carts. The priority is defined as the
    completion of the plate plugged in the cart times the priority its
    priority."""

    from sdss.internal.manga.Totoro.dbclasses import Plate

    session = db.Session()
    with session.begin():
        activePluggings = session.query(db.plateDB.ActivePlugging).all()

    priorities = []

    for cart in carts:
        cartActivePluggings = [aP for aP in activePluggings
                               if aP.plugging.cartridge.number == cart]
        if len(cartActivePluggings) == 0:
            raise exceptions.TotoroError('no active plugging for cart #{0}'
                                         .format(cart))

        priority = (cartActivePluggings[0].plugging
                    .plate.plate_pointings[0].priority)
        if priority == 0:
            priorities.append(0)
            continue

        try:
            completion = Plate(cartActivePluggings[0].plugging.plate.plate_id,
                               format='plate_id').getPlateCompletion()
        except exceptions.TotoroError:
            # If plate is not MaNGA (i.e., is APOGEE using one of our carts)
            completion = 0.  # Forces priority to be zero.

        priorities.append(completion * priority)

    return np.array(priorities)


def getCartPlate(cartNumber):
    """Returns the plate plugged in a cart or None."""

    session = db.Session()
    with session.begin():
        activePluggings = session.query(db.plateDB.ActivePlugging).all()

    for aP in activePluggings:
        if aP.plugging.cartridge.number == cartNumber:
            return aP.plugging.plate
    return None


class PluggerScheduler(object):
    """A class to schedule plugging requests."""

    def __init__(self, plates, jd0, jd1, **kwargs):

        log.debug('creating PluggerScheduler instance with JD0={0:.5f}, '
                  'JD1={1:.5f}'.format(jd0, jd1))

        self.timeline = Timeline(jd0, jd1, **kwargs)

        # Rejects plates with priority 1
        self._platesAtAPO = [plate for plate in plates
                             if plate.priority >
                             config['plugger']['noPlugPriority']]

        self.carts = OrderedDict([(key, None) for key in config['carts']])

        self._scheduleForced()

        if (len(self.timeline.plates) >= len(self.carts) or
                self.timeline.remainingTime <= 0):
            pass
        else:
            self.timeline.schedule([plate for plate in self._platesAtAPO
                                    if plate not in self.timeline.plates])

        self.allocateCarts(plates=self.timeline.plates)

        remainingTime = self.timeline.remainingTime
        if remainingTime > 0:
            log.important('{0:.2h} hours not allocated'.format(remainingTime))
        else:
            log.debug('all time has been allocated.')

    def _scheduleForced(self):
        """Schedules plates that have priority=10."""

        forcePlugPriority = int(config['plugger']['forcePlugPriority'])
        forcePlugPlates = [plate for plate in self._platesAtAPO
                           if plate.priority == forcePlugPriority]

        if len(forcePlugPlates) == 0:
            return

        # Sort forced plates so that the plates that are already plugged are
        # first.
        forcePlugPlatesSorted = [plate for plate in forcePlugPlates
                                 if plate.isPlugged] + \
                                [plate for plate in forcePlugPlates
                                 if not plate.isPlugged]

        if len(forcePlugPlatesSorted) > len(self.carts):
            warnings.warn('{0} plates marked with priority {1:d} but only '
                          '{2} carts available. Using the first {2} plates.'
                          .format(len(forcePlugPlatesSorted),
                                  forcePlugPriority, len(self.carts)),
                          exceptions.TotoroPluggerWarning)
            forcePlugPlatesSorted = forcePlugPlatesSorted[0:len(self.carts)]

        self.timeline.schedule(plates=forcePlugPlatesSorted, mode='plugger',
                               force=True)

        # Manually adds all the forced plates to the timeline, whether they
        # have been observed or not.
        self.timeline.allocateJDs(forcePlugPlatesSorted)
        self.timeline.plates += [plate for plate in forcePlugPlatesSorted
                                 if plate not in self.timeline.plates]

    def allocateCarts(self, plates=None, **kwargs):

        def logCartAllocation(cartNumber, plate, messages=None):
            if isinstance(messages, basestring):
                messages = [messages]
            plateid = plate.plate_id
            status = plate.statuses[0].label
            if status == 'Shipped' and plate.location.label == 'APO':
                msg = 'plate has not been marked'
                if messages is None:
                    messages = [msg]
                else:
                    messages.append(msg)
            msgStr = '({0})'.format(', '.join(messages)) \
                if messages is not None else ''
            log.important('Cart #{0} -> plate_id={1} {2}'
                          .format(cartNumber, plateid, msgStr))

        plates = plates if plates is not None else self._platesToAllocate

        mjd = int(self.timeline.endDate - 2400000.5)
        log.important('Plugging allocation for MJD={0:d} follows.'
                      .format(mjd))

        if len(plates) > len(self.carts):
            warnings.warn('{0} plates to allocate but only {1} carts '
                          'available. Using the first {1} plates.'
                          .format(len(plates), len(self.carts)),
                          exceptions.TotoroPluggerWarning)
            plates = plates[0:len(self.carts)]

        allocatedPlates = []

        for plate in plates:

            if plate.isPlugged:

                cartNumber = plate.getActiveCartNumber()

                self.carts[cartNumber] = plate.plate_id

                logCartAllocation(cartNumber, plate, 'already plugged')

                allocatedPlates.append(plate)

            else:

                cartNumber, status = getAvailableCart(
                    [cartNumber for cartNumber in self.carts
                     if self.carts[cartNumber] is None])

                if cartNumber is not None:

                    self.carts[cartNumber] = plate.plate_id

                    if status == 'empty':
                        logCartAllocation(cartNumber, plate, 'empty cart')
                    elif status == 'noMaNGAPlate':
                        logCartAllocation(cartNumber, plate,
                                          'replacing no-MaNGA plate')
                    elif status == 'complete':
                        logCartAllocation(cartNumber, plate,
                                          'replacing complete plate')
                    elif status == 'notStarted':
                        logCartAllocation(cartNumber, plate,
                                          'replacing non stated plate')

                    allocatedPlates.append(plate)

        unallocatedPlates = [plate for plate in plates
                             if plate not in allocatedPlates]
        unallocatedCarts = [cart for cart in self.carts
                            if self.carts[cart] is None]

        cartPriority = getCartPriority(unallocatedCarts)

        for plate in unallocatedPlates:
            bestCartIdx = np.argmin(cartPriority)
            cart = unallocatedCarts[bestCartIdx]
            self.carts[cart] = plate.plate_id
            logCartAllocation(cart, plate, 'replacing MaNGA plate')

            cartPriority = np.delete(cartPriority, bestCartIdx)
            unallocatedCarts.pop(bestCartIdx)

        if len(unallocatedCarts) > 0:
            for cart in unallocatedCarts:
                plate = getCartPlate(cart)
                if plate is None:
                    log.important('Cart #{0} -> no plate'.format(cart))
                else:
                    self.carts[cart] = plate.plate_id
                    logCartAllocation(cart, plate, 'kept plugged')
