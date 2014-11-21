#!/usr/bin/env python
# encoding: utf-8
"""
pluger.py

Created by José Sánchez-Gallego on 20 Oct 2014.
Licensed under a 3-clause BSD license.

Revision history:
    20 Oct 2014 J. Sánchez-Gallego
      Initial version
    15 Nov 2014 J. Sánchez-Gallego
      Improved the logic and added some convenience functions

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
    with session.begin(subtransactions=True):
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


def getCartStatus(cartNumber):
    """Returns the status of a cart."""

    return getAvailableCart([cartNumber])[1]


def getCartPriority(carts):
    """Returns the priority of a list of carts. The priority is defined as the
    completion of the plate plugged in the cart times the priority its
    priority."""

    from sdss.internal.manga.Totoro.dbclasses import Plate

    session = db.Session()
    with session.begin(subtransactions=True):
        activePluggings = session.query(db.plateDB.ActivePlugging).all()

    priorities = []

    for cart in carts:
        cartActivePluggings = [aP for aP in activePluggings
                               if aP.plugging.cartridge.number == cart]
        if len(cartActivePluggings) == 0:
            log.debug('no active plugging for cart #{0}'.format(cart))
            priorities.append(10)
            continue

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
    with session.begin(subtransactions=True):
        activePluggings = session.query(db.plateDB.ActivePlugging).all()

    for aP in activePluggings:
        if aP.plugging.cartridge.number == cartNumber:
            return aP.plugging.plate
    return None


def getCartForReplug(plate):
    """Returns the cart of the last plugging."""

    if len(plate.pluggings) == 0:
        return None

    scanMJDs = [plugging.fscan_mjd for plugging in plate.pluggings]
    return plate.pluggings[np.argmax(scanMJDs)].cartridge.number


class PluggerScheduler(object):
    """A class to schedule plugging requests."""

    def __init__(self, plates, jd0, jd1, **kwargs):

        log.debug('creating PluggerScheduler instance with JD0={0:.5f}, '
                  'JD1={1:.5f}'.format(jd0, jd1))

        self.timeline = Timeline(jd0, jd1, **kwargs)

        self.platesToSchedule = self.selectPlates(plates)
        log.info('scheduling {0} plates'.format(len(self.platesToSchedule)))

        self.carts = OrderedDict([(key, None) for key in config['carts']])

        self._scheduleForced()

        if (len(self.timeline.plates) >= len(self.carts) or
                self.timeline.remainingTime <= 0):
            pass
        else:
            self.timeline.schedule([plate for plate in self.platesToSchedule
                                    if plate not in self.timeline.plates],
                                   mode='plugger')

        self.allocateCarts(plates=self.timeline.plates)

        remainingTime = self.timeline.remainingTime
        if remainingTime > 0:
            log.important('{0:.2f}h hours not allocated'.format(remainingTime))
        else:
            log.debug('all time has been allocated.')

    def selectPlates(self, plates):
        """Selects plates to schedule, rejecting those which are invalid or
        must not be scheduled."""

        prioPlug = int(config['plugger']['forcePlugPriority'])

        # Selects plates with priority 10 or plugged
        platesToSchedule = []

        for plate in plates:
            if (plate.priority == prioPlug or plate.isPlugged):
                platesToSchedule.append(plate)

        # Adds the remainder of the plates
        for plate in [plate for plate in plates
                      if plate not in platesToSchedule]:

            plate_id = plate.plate_id

            if plate.priority <= config['plugger']['noPlugPriority']:
                log.debug('Skipped plate_id={0} because of low priority'
                          .format(plate_id))
                continue

            if plate.isComplete:
                log.debug('Skipped plate_id={0} because is complete'
                          .format(plate_id))
                continue

            # If the plate has been started but is not plugged and the cart
            # that must be used is already plugged, skips the plate.
            if plate.getPlateCompletion(includeIncompleteSets=True) > 0:
                cartToUse = getCartForReplug(plate)
                isCartFree = True
                for pp in platesToSchedule:
                    if (pp.isPlugged and
                            pp.getActiveCartNumber() == cartToUse):
                        isCartFree = False
                        break
                if not isCartFree:
                    log.info('Skipped plate_id={0} because is a replug and '
                             'its cart is in use.'.format(plate_id))
                    continue
                else:
                    # Marks the plate for future reference
                    plate.isReplug = True

            platesToSchedule.append(plate)

        return platesToSchedule

    def _scheduleForced(self):
        """Schedules plates that have priority=10."""

        forcePlugPriority = int(config['plugger']['forcePlugPriority'])
        forcePlugPlates = [plate for plate in self.platesToSchedule
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

        self.timeline.schedule(plates=forcePlugPlatesSorted, mode='plugger')

        # Manually adds all the forced plates to the timeline, whether they
        # have been observed or not.
        self.timeline.allocateJDs(forcePlugPlatesSorted)
        self.timeline.plates += [plate for plate in forcePlugPlatesSorted
                                 if plate not in self.timeline.plates]

    def allocateCarts(self, plates, **kwargs):

        def logCartAllocation(cartNumber, plate, messages=''):
            """Convenience function to log the cart allocation."""

            if plate is None:
                log.important('Cart #{0} -> empty'.format(cartNumber))
                return

            if not isinstance(messages, (list, tuple)):
                messageList = [messages]

            plateid = plate.plate_id
            status = plate.statuses[0].label

            if status == 'Shipped' and plate.location.label == 'APO':
                messageList.append('plate has not been marked')

            if hasattr(plate, 'isReplug') and plate.isReplug:
                messageList.append('replug')

            jointMessage = ', '.join(messageList)
            msgStr = '({0})'.format(jointMessage) if jointMessage != '' else ''

            log.important('Cart #{0} -> plate_id={1} {2}'
                          .format(cartNumber, plateid, msgStr))

        if len(plates) > len(self.carts):
            warnings.warn('{0} plates to allocate but only {1} carts '
                          'available. Using the first {1} plates.'
                          .format(len(plates), len(self.carts)),
                          exceptions.TotoroPluggerWarning)
            plates = plates[0:len(self.carts)]

        mjd = int(self.timeline.endDate - 2400000.5)
        log.important('Plugging allocation for MJD={0:d} follows:'
                      .format(mjd))

        allocatedPlates = []

        # Sorts plates so that plates to replug are the first ones to be
        # allocated
        sortedPlates = [plate for plate in plates
                        if hasattr(plate, 'isReplug') and plate.isReplug]
        sortedPlates += [plate for plate in plates
                         if plate not in sortedPlates]

        # Dictionary to save the cart allocation
        cartPlateMessage = {}

        for plate in plates:

            if plate.isPlugged:

                cartNumber = plate.getActiveCartNumber()
                self.carts[cartNumber] = plate.plate_id
                logCartAllocation(cartNumber, plate, 'already plugged')
                allocatedPlates.append(plate)

            else:

                if hasattr(plate, 'isReplug') and plate.isReplug:
                    cartNumber = getCartForReplug(plate)
                    status = getCartStatus(cartNumber)
                else:
                    cartNumber, status = getAvailableCart(
                        [cartNumber for cartNumber in self.carts
                         if self.carts[cartNumber] is None])

                if cartNumber is not None:

                    self.carts[cartNumber] = plate.plate_id

                    if status == 'empty':
                        cartPlateMessage[cartNumber] = (plate, 'empty cart')
                    elif status == 'noMaNGAPlate':
                        cartPlateMessage[cartNumber] = (
                            plate, 'replacing no-MaNGA plate')
                    elif status == 'complete':
                        cartPlateMessage[cartNumber] = (
                            plate, 'replacing complete plate')
                    elif status == 'notStarted':
                        cartPlateMessage[cartNumber] = (
                            plate, 'replacing non started plate')

                    allocatedPlates.append(plate)

        # Now it handles the plates that will require replacing a MaNGA plate
        unallocatedPlates = [plate for plate in plates
                             if plate not in allocatedPlates]
        unallocatedCarts = [cart for cart in self.carts
                            if self.carts[cart] is None]

        cartPriority = getCartPriority(unallocatedCarts)

        for plate in unallocatedPlates:
            bestCartIdx = np.argmin(cartPriority)
            cart = unallocatedCarts[bestCartIdx]
            self.carts[cart] = plate.plate_id
            cartPlateMessage[cart] = (plate, 'replacing MaNGA plate')

            cartPriority = np.delete(cartPriority, bestCartIdx)
            unallocatedCarts.pop(bestCartIdx)

        # If there are unallocated carts, leaves them untouched.
        if len(unallocatedCarts) > 0:
            for cart in unallocatedCarts:
                plate = getCartPlate(cart)
                if plate is None:
                    cartPlateMessage[cart] = (None, '')
                else:
                    self.carts[cart] = plate.plate_id
                    cartPlateMessage[cart] = (plate, 'unchanged')

        # Logs the allocation
        for cart in sorted(cartPlateMessage.keys()):
            logCartAllocation(cart, cartPlateMessage[cart][0],
                              cartPlateMessage[cart][1])
