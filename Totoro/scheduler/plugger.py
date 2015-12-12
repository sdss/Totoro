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
from Totoro import log, config, site
from Totoro.db import getConnection
from Totoro.scheduler.timeline import Timeline
from Totoro.exceptions import TotoroPluggerWarning, TotoroPluggerError
from Totoro.utils import intervals
from collections import OrderedDict
import warnings
import numpy as np


__all__ = ['Plugger']

cartStatusCodes = {0: 'empty', 1: 'noMaNGAplate', 2: 'MaNGA_complete',
                   3: 'MaNGA_noStarted', 4: 'MaNGA_started', 10: 'unknown'}

replaceMsgs = {0: 'empty cart', 1: 'replacing non-MaNGA plate',
               2: 'replacing complete MaNGA plate',
               3: 'replacing non-started MaNGA plate',
               4: 'replacing started MaNGA plate',
               10: 'replacing plate with unknown status'}


def getForcePlugPlates():
    """Returns a list of plates with priority `forcePlugPriority`."""

    from Totoro.dbclasses.plate import Plates

    forcePlugPriority = int(config['plugger']['forcePlugPriority'])

    db = getConnection()
    session = db.Session()
    plateDB = db.plateDB

    with session.begin():
        plates = session.query(plateDB.Plate).join(
            plateDB.PlateToSurvey, plateDB.Survey, plateDB.SurveyMode,
            plateDB.PlatePointing, plateDB.PlateLocation).filter(
                plateDB.Survey.label == 'MaNGA',
                plateDB.SurveyMode.label == 'MaNGA dither',
                plateDB.PlateLocation.label == 'APO',
                plateDB.PlatePointing.priority == forcePlugPriority).order_by(
                    plateDB.Plate.plate_id).all()

    return Plates(plates)


def getActivePluggings():
    """Returns a list with the active pluggings."""

    db = getConnection()
    session = db.Session()

    # Gets active pluggings
    with session.begin():
        activePluggings = session.query(
            db.plateDB.ActivePlugging).order_by(
                db.plateDB.ActivePlugging.pk).all()

    return activePluggings


def getCartStatus(activePluggings, cartNumber):
    """Returns the status of the plate in a cart.

    The returned tuple is (cart_number, plate, status_code, completion)
    """

    from Totoro.dbclasses.plate import Plate

    cartActivePluggings = [aP for aP in activePluggings
                           if aP.plugging.cartridge.number == cartNumber]

    if len(cartActivePluggings) == 0:
        return (cartNumber, None, 0, 0)  # Empty cart
    elif len(cartActivePluggings) > 1:
        raise TotoroPluggerError(
            'PLUGGER: something went wrong. Cart #{0} has more than one '
            'active plugging'.format(cartNumber))
    else:
        cartActivePlugging = cartActivePluggings[0]

    plate = cartActivePlugging.plugging.plate

    isMaNGAPlate = (plate.currentSurveyMode is not None and
                    'MaNGA' in plate.currentSurveyMode.label)

    if not isMaNGAPlate:
        return (cartNumber, plate, 1, 0)  # No MaNGA plate

    totoroPlate = Plate(plate)

    if totoroPlate.isComplete:
        return (cartNumber, totoroPlate, 2, 1.)  # Complete MaNGA plate

    plateCompletion = totoroPlate.getPlateCompletion()
    if plateCompletion == 0:
        # Non-stated MaNGA plate
        return (cartNumber, totoroPlate, 3, plateCompletion)
    else:
        # Started MaNGA plate
        return (cartNumber, totoroPlate, 4, plateCompletion)

    return (cartNumber, None, 10, 0)  # Unknown


def getCartPlate(activePluggings, cartNumber):
    """Returns the plate plugged in a cart or None."""

    for aP in activePluggings:
        if aP.plugging.cartridge.number == cartNumber:
            return aP.plugging.plate
    return None


def getCartLastPlugging(plate):
    """Returns the cart number of the last plugging."""

    if len(plate.pluggings) == 0:
        return None

    scanMJDs = [plugging.fscan_mjd for plugging in plate.pluggings]
    return plate.pluggings[np.argmax(scanMJDs)].cartridge.number


def prioritiseCarts(carts):
    """Returns a list of carts sorted by priority for being allocated.

    Parameters
    ----------
    carts : list
        A list with all or a subset of the MaNGA carts. Each element in the
        list is a tuple of the form
        `(cart_number, plate, status_code, completion)`

    Returns
    -------
    result : list
        Returns `carts` sorted in the following order:
        - Empty carts
        - Carts with complete plates
        - Carts with plates of unknown status
        - Carts with non-MaNGA plates
        - Carts with MaNGA plates that have not been started
        - Carts with started plates sorted by completion.

    """

    # Creates some intermediate lists
    empty = []
    noMaNGA = []
    complete = []
    noStarted = []
    unknown = []
    started = []

    # Assigns carts to the appropriate list.
    for cart in carts:
        statusLabel = cartStatusCodes[cart[2]]
        if statusLabel == 'empty':
            empty.append(cart)
        elif statusLabel == 'noMaNGAplate':
            noMaNGA.append(cart)
        elif statusLabel == 'MaNGA_complete':
            complete.append(cart)
        elif statusLabel == 'MaNGA_noStarted':
            noStarted.append(cart)
        elif statusLabel == 'unknown':
            unknown.append(cart)
        elif statusLabel == 'MaNGA_started':
            started.append(cart)

    # Sorts started carts using the completion (fourth element of the tuple).
    started = sorted(started, key=lambda xx: xx[3])

    # Returns carts in the desired order.
    return empty + complete + unknown + noMaNGA + noStarted + started


class Plugger(object):
    """A class to schedule plugging requests.

    A new `Plugger` instance is initiated by providing a `startDate` and
    `endDate`. If both of them are None, no scheduling is performed and only
    the already plugged plates are included.

    Parameters
    ----------
    startDate : float or None
        The JD at which the scheduling begins.
    endDate : float or None
        The JD at which the scheduling ends.

    """

    def __init__(self, startDate=None, endDate=None, **kwargs):

        startDate = None if not startDate else startDate
        endDate = None if not endDate else endDate

        self._nNewExposures = []
        self._platesToSchedule = []

        # Runs init method depending on startDate and endDate.
        if startDate is None and endDate is None:
            self._initNoManga()
        elif any([startDate, endDate]) and not all([startDate, endDate]):
            raise TotoroPluggerError('PLUGGER: either startDate = endDate = '
                                     'None or both need to be defined.')
        else:
            self._initFromDates(startDate, endDate)

    def _initNoManga(self):
        """Inits a Plugger instance when no MaNGA time is scheduled.

        In this case, the cart assignement contains only those plugged MaNGA
        plates that are not complete, sorted by preference of the carts being
        overridden by APOGEE."""

        from Totoro.dbclasses import getPlugged

        self.startDate = None
        self.endDate = None

        warnings.warn('PLUGGER: no JD1, JD2 values provided. Plugger will '
                      'only return plugged, on-completed plates.',
                      TotoroPluggerWarning)

        # Get MaNGA plugged plates
        pluggedPlates = getPlugged(fullCheck=False, updateSets=False)

        self.carts = OrderedDict()
        self._nNewExposures = dict()

        for plate in pluggedPlates:
            if not plate.isComplete:
                cart = plate.getActiveCartNumber()
                self.carts[cart] = plate

    def _initFromDates(self, jd0, jd1, **kwargs):
        """Initialites the Plugger instance from two JD dates.

        This method does not actually schedules plates for the range
        `[jd0, jd1]`. Instead, it creates the `Totoro.Timeline` object for
        this plugging requests and obtains the list of plates that can be
        scheduled. The real scheduling happens when `Plugger.schedule()`
        is called.

        """

        assert jd0 < jd1, 'JD1 cannot be lower than JD1'

        self.startDate = jd0
        self.endDate = jd1
        scheduledTime = (self.endDate - self.startDate) * 24.

        log.info('PLUGGER: start date: {0}'.format(self.startDate))
        log.info('PLUGGER: end date: {0}'.format(self.endDate))
        log.info('PLUGGER: scheduling {0:.2f} hours'.format(scheduledTime))

        # Creates the timeline object for this plugging request.
        self.timeline = Timeline(self.startDate, self.endDate)

        # Determines the plates to schedule.
        self._platesToSchedule = self.getPlatesToSchedule(**kwargs)
        log.info('PLUGGER: scheduling {0} plates'
                 .format(len(self._platesToSchedule)))

        # Initialises a dictionary with the MaNGA carts.
        self.carts = OrderedDict([(key, None) for key in config['mangaCarts']])

    def getPlatesToSchedule(
            self, onlyMarked=False,
            onlyVisiblePlates=config['plugger']['onlyVisiblePlates']):
        """Selects plates to schedule.

        Determines the list of plates to schedule by rejecting those which
        are invalid or outside the LST window for the night.

        """

        from Totoro import dbclasses

        assert isinstance(onlyVisiblePlates, int), \
            'onlyVisiblePlates must be a boolean'

        log.info('PLUGGER: getting plates at APO with onlyMarked={0}'
                 .format(onlyMarked))

        # If we are only selecting plates observable that night, determines
        # the RA range of the plates to accept.
        if onlyVisiblePlates:
            lstRange = site.localSiderealTime([self.startDate, self.endDate])
            window = config['plateVisibilityMaxHalfWindowHours']
            raRange = np.array([(lstRange[0] - window) * 15.,
                                (lstRange[1] + window) * 15.])

            log.info('PLUGGER: selecting plates with RA in range {0}'
                     .format(str(raRange % 360)))

            # If the RA range wraps around 0, we split it in two non-wrapping
            # ranges
            raRange = intervals.splitInterval(raRange, 360.)

        else:
            raRange = None

        # Selects plates at APO with the appropriate parameters.
        platesAtAPO = dbclasses.getAtAPO(onlyIncomplete=True,
                                         onlyMarked=onlyMarked,
                                         rejectLowPriority=True,
                                         fullCheck=False,
                                         updateSets=False,
                                         raRange=raRange)
        plugged = dbclasses.getPlugged()

        platesToSchedule = platesAtAPO + [plate for plate in plugged
                                          if plate not in platesAtAPO]

        log.info('PLUGGER: plates found: {0}'.format(len(platesToSchedule)))

        return platesToSchedule

    def schedule(self, **kwargs):
        """Schedules the selected plates and gets the list of plates to plug.

        Schedules the plates selected during the initialisation of the instance
        and determines the list of of plates to plug and the cart allocation.

        """

        forcePlugPriority = int(config['plugger']['forcePlugPriority'])

        # Removes force-plug plates from the list of plates to schedule.
        # We'll add them back at the end, but we don't use them to cover the
        # scheduled time.
        self._platesToSchedule = [plate for plate in self._platesToSchedule
                                  if plate.priority < forcePlugPriority]

        # Gets a list of force-plug plates
        forcePlugPlates = getForcePlugPlates()

        # If there are more force plug plates than carts, there is no point in
        # scheduling the rest.
        if len(forcePlugPlates) < len(self.carts):
            self.timeline.schedule(self._platesToSchedule, mode='plugger',
                                   **kwargs)

        # Now we add back the force plug plates
        scheduledPlates = self.timeline.scheduled + forcePlugPlates

        # We log the number of new exposures for the plates in the timeline.
        # We'll use this later when we prioritise carts.
        self._nNewExposures = {}
        for plate in scheduledPlates:
            nNewExposuresPlate = len(plate.getMockExposures())
            self._nNewExposures[plate.plate_id] = nNewExposuresPlate

        # Allocates carts
        self.allocateCarts(scheduledPlates)

        remainingTime = self.timeline.remainingTime
        if remainingTime > 0:
            log.important('PLUGGER: {0:.2f}h hours not allocated'
                          .format(remainingTime))
        else:
            log.debug('PLUGGER: all the time has been allocated.')

    def logCartAllocation(self, activePluggings):
        """Logs the cart allocation."""

        mjd = int(self.timeline.endDate - 2400000.5)
        log.important('Plugging allocation for MJD={0:d} follows:'
                      .format(mjd))

        for cartNo, plate in self.carts.iteritems():

            cartStatus = getCartStatus(activePluggings, cartNo)
            pluggedPlate = cartStatus[1]
            status = cartStatusCodes[cartStatus[2]]

            if plate is None:

                if status == 'noMaNGAplate':
                    message = ('plate_id={0} (APOGEE-2 plate, not doing '
                               'anything)'.format(pluggedPlate.plate_id))
                elif cartNo in config['offlineCarts']:
                    message = 'offline'
                elif pluggedPlate is not None:
                    message = ('plate_id={0} (unplug)'
                               .format(pluggedPlate.plate_id))
                else:
                    message = 'empty'

                log.important('Cart #{0} -> {1}'
                              .format(cartNo, message))

            else:

                if (pluggedPlate is not None and
                        pluggedPlate.plate_id == plate.plate_id):
                    message = 'already plugged'
                else:
                    message = replaceMsgs[cartStatus[2]]

                    if getCartLastPlugging(plate) is not None:
                        message += ', replug'

                    plateStatus = plate.statuses[0].label

                    if (plateStatus == 'Shipped' and
                            plate.location.label == 'APO'):
                        message += ', plate has not been marked'

                log.important('Cart #{0} -> plate_id={1} ({2})'
                              .format(cartNo, plate.plate_id, message))

    def _getCart(self, sortedCarts):
        """Given a list of sorted carts returns the first not allocated."""

        for cart in sortedCarts:
            if self.carts[cart[0]] is None:
                return cart

    def allocateCarts(self, plates):
        """Allocates plates into carts in the most efficient way."""

        offlineCarts = config['offlineCarts']
        forcePlugPriority = int(config['plugger']['forcePlugPriority'])

        if len(plates) > len(self.carts):
            warnings.warn('PLUGGER: {0} plates to allocate but only {1} carts '
                          'available. Using the first {1} plates.'
                          .format(len(plates), len(self.carts)),
                          TotoroPluggerWarning)
            plates = plates[0:len(self.carts)]

        # Gets the active pluggings
        activePluggings = getActivePluggings()

        # Gets the status of the plates in each cart.
        cartStatus = [getCartStatus(activePluggings, cartNumber)
                      for cartNumber in self.carts
                      if cartNumber not in offlineCarts]

        # Sorts carts by priority.
        sortedCarts = prioritiseCarts(cartStatus)

        # Gets a list of plugged plates and the carts in which they are plugged
        plugged = [plate for plate in plates if plate.isPlugged]
        pluggedCarts = [plate.getActiveCartNumber() for plate in plugged]

        allocatedPlates = []

        # Allocates force-plug plates. If the plate has been plugged before
        # tries to use the same cart, unless that cart is offline or contains
        # a plate that we want to keep plugged.
        forcePlugPlates = [plate for plate in plates
                           if plate.priority == forcePlugPriority]

        for plate in forcePlugPlates:
            lastCart = getCartLastPlugging(plate)
            if (lastCart is not None and lastCart not in offlineCarts and
                    lastCart not in pluggedCarts):
                self.carts[lastCart] = plate
            else:
                cartData = self._getCart(sortedCarts)
                self.carts[cartData[0]] = plate
            allocatedPlates.append(plate)

        # Allocates plates that are already plugged.
        for plate in plugged:
            if plate in allocatedPlates:
                continue
            cartNumber = plate.getActiveCartNumber()
            if self.carts[cartNumber] is None:
                self.carts[cartNumber] = plate
            else:
                cartData = self._getCart(sortedCarts)
                self.carts[cartData[0]] = plate
            allocatedPlates.append(plate)

        # Allocates replugs.
        for plate in plates:
            if plate in allocatedPlates:
                continue

            lastCart = getCartLastPlugging(plate)
            if lastCart is None or lastCart in offlineCarts:
                continue

            if self.carts[lastCart] is None:
                self.carts[lastCart] = plate
                allocatedPlates.append(plate)
            else:
                log.debug('PLUGGER: not plugging plate {0} in its '
                          'original cart {1} because it is not available'
                          .format(plate.plate_id, lastCart))
                continue

        # Allocates the remaining plates
        for plate in plates:
            if plate in allocatedPlates:
                continue

            cartData = self._getCart(sortedCarts)
            cartNumber, pluggedPlate, statusCode, completion = cartData
            self.carts[cartNumber] = plate
            allocatedPlates.append(plate)

        if len(plates) > len(allocatedPlates):
            warnings.warn('PLUGGER: {0} plates have not been allocated'
                          .format(len(plates) - len(allocatedPlates)),
                          TotoroPluggerWarning)

        remainingCarts = [cart for cart in sortedCarts
                          if self.carts[cart[0]] is None]

        # Checks unassigned carts
        for cart in remainingCarts:
            cartNumber, pluggedPlate, statusCode, completion = cart
            if completion >= 1:
                continue
            elif (cartStatusCodes[statusCode] != 'noMaNGAplate' and
                    pluggedPlate is not None):
                # If this is a MaNGA plate, keeps it.
                self.carts[cartNumber] = pluggedPlate

        # Logs the allocation
        self.logCartAllocation(activePluggings)

    def getASOutput(self, **kwargs):
        """Returns the plugging request in the autoscheduler format."""

        if self.startDate is not None and self.endDate is not None:
            mode = 'mangaLead'
            self.schedule(**kwargs)
        else:
            mode = 'apogeeLead'

        # Now we add a list with the priority order of the allocated carts.
        # This is useful if APOGEE needs to take over some of our carts.
        # In this way, they'll first use our carts with lower priority.
        cartOrder = self.getCartOrder(mode=mode)

        # Removes cart without an allocated MaNGA plate
        carts = OrderedDict([(key, value) for key, value in self.carts.items()
                             if value is not None])

        # First we add carts not used to cart_order, with lower priority
        nonUsedCarts = [cartNo for cartNo in config['mangaCarts']
                        if cartNo not in cartOrder]

        cartOrder = nonUsedCarts + cartOrder

        # We also add the APOGEE carts
        cartOrder = config['apogeeCarts'][::-1] + cartOrder

        # We change the Totoro.Plate instances to plate_ids
        for key in carts:
            if key != 'cart_order' and carts[key] is not None:
                carts[key] = carts[key].plate_id

        carts['cart_order'] = cartOrder

        return carts

    def getCartOrder(self, mode='mangaLead'):
        """Returns a prioritised cart list.

        If `mode='mangaLead'`, priority (from low to high) goes as it follows:

        (1) Completed plates. Note that there should be no completed plates
            in self.carts, but this is just to double-check.
        (2) Other plates, sorted by the number of exposures needed to fill out
            the night.
        (3) Force plug plates (plates with priority 10)

        If `mode='apogeeLead'`, (2) is replaced with the number of incomplete
        sets. Plates with incomplete sets are given the highest priority. This
        is used when `addCartOrder` is called for a `Plugger` object
        initialised without dates, and not scheduling is performed.

        Additionally, we give high priority to offline cart. This is to avoid
        APOGEE to plug co-designed plates in those carts if possible. If
        `mode='mangaLead'``, offline carts are given the highest priority after
        all scheduled carts. If ``mode='apogeeLead'`` offline carts are given
        higher priority than any but carts with incomplete sets.

        """

        forcePlugPriority = int(config['plugger']['forcePlugPriority'])

        assert mode in ['apogeeLead', 'mangaLead'], \
            'metric must be "apogeeLead" or "mangaLead"'

        # Splits plates in the three priority categories. We use useMock=False
        # because these plates contain mock exposures from the simulation.

        completed = []
        scheduled = []
        forcePlug = []

        for cart, plate in self.carts.iteritems():
            if plate is None:
                continue
            if plate.priority == forcePlugPriority:
                forcePlug.append((cart, plate))
            elif plate.getPlateCompletion(useMock=False) > 1:
                completed.append((cart, plate))
            else:
                scheduled.append((cart, plate))

        if mode == 'mangaLead':

            # Retrieves how many scheduled (mock) exposures are in each plate.
            nExposures = [self._nNewExposures[plate.plate_id]
                          if plate.plate_id in self._nNewExposures else 0
                          for cart, plate in scheduled]

            # Sorts scheduled exposures from few to many scheduled exposures.
            scheduledOrdered = [scheduled[ii] for ii in np.argsort(nExposures)]

        elif mode == 'apogeeLead':

            # Finds out what plates have incomplete sets.
            incompleteSets = []
            completeSets = []

            for ii in range(len(scheduled)):
                if scheduled[ii][1].hasIncompleteSets():
                    incompleteSets.append(scheduled[ii])
                else:
                    completeSets.append(scheduled[ii])

            # Sorts the scheduled plates according to completion. For plates
            # with incomplete sets we take them into account.

            sortedIncompleteSets = sorted(
                incompleteSets, key=lambda xx: xx[1].getPlateCompletion(
                    includeIncompleteSets=True))

            sortedCompleteSets = sorted(
                completeSets, key=lambda xx: xx[1].getPlateCompletion())

            # We put plates with complete sets first
            scheduledOrdered = sortedCompleteSets + sortedIncompleteSets

        usedCarts = [cart for cart, plate in
                     completed + scheduledOrdered + forcePlug]

        # Creates master ordered list
        if mode == 'mangaLead':
            offline = [(cart, None) for cart in config['offlineCarts']
                       if cart not in usedCarts]
            orderedCarts = completed + offline + scheduledOrdered + forcePlug
        else:
            # Identifies the first plate with incomplete sets
            ii = 0
            for cart, plate in scheduledOrdered:
                if np.any([ss.getStatus()[0] == 'Incomplete'
                           for ss in plate.sets]):
                    break
                ii += 1

            # Adds offline carts before
            for cart in config['offlineCarts']:
                if cart in usedCarts:
                    continue
                scheduledOrdered.insert(ii, (cart, None))

            orderedCarts = completed + scheduledOrdered + forcePlug

        return [cart for cart, plate in orderedCarts]
