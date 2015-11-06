# Change Log

## [1.0.1] - 2015-11-06
### Added
- New method in Exposure to get the mean airmass during the exposure.
### Changed
- Zenith avoidance set to 5 degrees by default for planner.
### Fixed
- Small bug fixes
- Fixes a bug when there are no exposures to clean in a simulation.

## [1.0.0] - 2015-10-30
### Added
- CHANGELOG.md
- SDSSconnect added to Totoro root
- EUPS
- Plugger: offline carts are now given higher cart_order. This is because we
don't want APOGEE to use MaNGA offline carts for co-designed plates, if
possible. When the plugger is running for a MaNGA night offline carts are
given the highest priority after all scheduled carts. For non-MaNGA nights,
offline carts are given higher priority than any other cart except carts with
plates with incomplete sets.
- Module file
### Changed
- Totoro moved to its own repo.
- scheduler_utils: plates are sorted by completion so that, if several plates
complete with the same number of exposures, the one with highest completion
gets chosen.
- Now uses SDSSconnect for connecting to the DB. Requires passwords to be
defined in .pgpass.
- DB classes now don't inherit from model classes. Instead, the DB query object
is recorded as object.\_dbObject. DB attributes can still be accessed from the
Totoro object via a custom `__getattr__` method.
