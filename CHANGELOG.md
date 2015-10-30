# Change Log

## [Unreleased]
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
- Now uses SDSSconnect for connecting to the DB. Requires passwords to be
defined in .pgpass.
- DB classes now don't inherit from model classes. Instead, the DB query object
is recorded as object.\_dbObject. DB attributes can still be accessed from the
Totoro object via a custom `__getattr__` method.
