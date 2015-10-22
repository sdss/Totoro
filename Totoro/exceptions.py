# Exceptions for pyPhotom classes and functions

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


class TotoroError(Exception):
    """Base exception for Totoro. Other exceptions should inherit this."""
    pass


class TotoroNotImplemented(TotoroError):
    """A class for exceptions about functionalities not yet implemented."""
    pass


class PlateNotPlugged(TotoroError):
    """Exception for when a plate is not currently plugged."""
    pass


class EmptySet(TotoroError):
    """Exception for when a set has no exposures."""
    pass


class TotoroPluggerError(TotoroError):
    """Exception for the Plugger class."""
    pass


class NoMangaPlate(TotoroError):
    """Exception to be raised if a plate is not MaNGA."""
    pass


class TotoroSubtransactionError(TotoroError):
    """Exception to be raised if Totoro is called from within a DB session."""
    pass


class TotoroWarning(Warning):
    """Base warning for Totoro."""
    pass


class TotoroUserWarning(UserWarning, TotoroWarning):
    """The primary warning class."""
    pass


class TotoroPluggerWarning(TotoroUserWarning):
    """Warning for Plugger issues."""
    pass


class TotoroPlannerWarning(TotoroUserWarning):
    """Warning for Planner issues."""
    pass


class DustMapWarning(TotoroUserWarning):
    """A warning for when no dust map is present."""
    pass


class FieldWarning(TotoroUserWarning):
    """A warning for Field."""
    pass


class MultiplePlatePointings(TotoroUserWarning):
    """Warning for when a plate has more than one plate pointing."""
    pass


class NoMangaExposure(TotoroUserWarning):
    """Warning in case an exposure has not mangaDB counterpart."""
    pass


class ExposureWithoutSet(TotoroUserWarning):
    """Warning in case an exposure has no set."""
    pass


class NoObservingBlock(TotoroUserWarning):
    """Warning in case a MJD has no MaNGA observing time."""
    pass
