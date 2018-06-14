# Totoro

This directory contains the MaNGA autoscheduler logic (codename Totoro) as well as classes
and functions for night planning and plugging.

This code has been created and is maintained by José Sánchez-Gallego (gallegoj@uw.edu)

## Dependences

Totoro requires Python 2.7+ (including Python 3) to work (previous versions of Python may work but the code has not been tested). The support for Python 3 is still tentative. If you find a bug please create an [issue](https://github.com/sdss/Totoro/issues/new).

Totoro depends on the following third-party libraries:

- astropy >= 2.0
- sqlalchemy
- pyyaml
- pydl

The following SDSS dependences are required:

- autoscheduler (if not present, a local copy of the schedule file will be used)

## License

Totoro is licensed under a 3-clause BSD style license.
