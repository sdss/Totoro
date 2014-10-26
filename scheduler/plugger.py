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
from collections import OrderedDict
import warnings
import numpy as np


db = TotoroDBConnection()


def getAvailableCart(carts):
    """Gets a cart that without a plugged plate. If none is available, returns
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

    cartStatus = OrderedDict([('complete', []),
                              ('noMaNGAPlate', []),
                              ('notStarted', [])])

    for aP in activePluggings:

        for key in cartStatus:
            if len(cartStatus[key]) > 0:
                return (cartStatus[key][0], key)

        cartNumber = aP.plugging.cartridge.number

        if cartNumber not in carts:
            continue

        if aP.plugging.plate.isComplete:
            cartStatus['complete'].append(cartNumber)

        elif (aP.plugging.plate.currentSurveyMode is None or
                'MaNGA' not in aP.plugging.plate.currentSurveyMode.label):
            cartStatus['noMaNGAPlate'].append(cartNumber)

        elif (Plate(aP.plugging.plate.plate_id, format='plate_id')
              .getPlateCompletion() == 0):
            cartStatus['notStarted'].append(cartNumber)

    return None


def getCartPriority(carts):
    """Returns the priority of a list of plates. The priority is defined as the
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
    """A class to schedule plugging recuests."""

    def __init__(self, plates, jd0, jd1, **kwargs):

        log.debug('creating PluggerScheduler instance with JD0={0:.5f}, '
                  'JD1={1:.5f}'.format(jd0, jd1))

        self.timeline = Timeline(jd0, jd1, **kwargs)
        self._platesAtAPO = plates
        self._platesToAllocate = []

        self.carts = OrderedDict([(key, None) for key in config['carts']])

        self._scheduleForced()

        # if (len(self._platesToAllocate) >= len(self.carts) or
        #         self.timeline.lstRange.size == 0):
        self.allocateCarts(plates=self._platesToAllocate)
        # else:
        #     self.timeline.schedule([plate for plate in self._platesAtAPO
        #                             if plate not in self._platesToAllocate])

    def _scheduleForced(self):

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

        print(forcePlugPlatesSorted[0].getPlateCompletion())
        self.timeline.schedule(forcePlugPlatesSorted, mode='plugger',
                               force=True)
        print(forcePlugPlatesSorted[0].sets)
        print([len(set.totoroExposures) for set in forcePlugPlatesSorted[0].sets])

        self._platesToAllocate += forcePlugPlates

    def allocateCarts(self, plates=None, **kwargs):

        plates = plates if plates is not None else self._platesToAllocate

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

                log.important('Cart #{0} -> plate_id={1} (already plugged)'
                              .format(cartNumber, plate.plate_id))

                allocatedPlates.append(plate)

            else:

                cartNumber, status = getAvailableCart(
                    [cartNumber for cartNumber in self.carts
                     if self.carts[cartNumber] is None])

                if cartNumber is not None:

                    self.carts[cartNumber] = plate.plate_id

                    if status == 'empty':
                        log.important('Cart #{0} -> plate_id={1} (empty cart)'
                                      .format(cartNumber, plate.plate_id))
                    elif status == 'noMaNGAPlate':
                        log.important('Cart #{0} -> plate_id={1} (replacing '
                                      'no-MaNGA plate)'.format(
                                          cartNumber, plate.plate_id))
                    elif status == 'complete':
                        log.important('Cart #{0} -> plate_id={1} (replacing '
                                      'complete plate)'.format(
                                          cartNumber, plate.plate_id))
                    elif status == 'notStarted':
                        log.important('Cart #{0} -> plate_id={1} (replacing '
                                      'not started plate)'.format(
                                          cartNumber, plate.plate_id))

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
            log.important('Cart #{0} -> plate_id={1} (replacing '
                          'MaNGA plate)'.format(cart, plate.plate_id))

            cartPriority = np.delete(cartPriority, bestCartIdx)
            unallocatedCarts.pop(bestCartIdx)

        if len(unallocatedCarts) > 0:
            for cart in unallocatedCarts:
                plate = getCartPlate(cart)
                if plate is None:
                    log.important('Cart #{0} -> no plate'.format(cart))
                else:
                    self.carts[cart] = plate.plate_id
                    log.important('Cart #{0} -> plate_id={1} (kept plugged)'
                                  .format(cart, plate.plate_id))
