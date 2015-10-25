# Totoro

This directory contains the MaNGA autoscheduler logic (codename Totoro) as well as classes
and functions for night planning and plugging.

## Author

This code has been created and is maintained by José Sánchez-Gallego (j.sanchezgallego@uky.edu)

## Dependences

Totoro requires Python 2.7 to work (previous versions of Python may work but the code has not been tested).
Python 3 can be used but is not recommended yet as not all the dependences have yet been updated.

Totoro depends on the following third-party libraries:

- astropy >= 1.0
- astropysics
- sqlalchemy
- pyyaml

The following SDSS dependences are required:

- sdss_python_module
- autoscheduler (if not present, a local copy of the schedule file will be used)

## License

Totoro is licensed under a 3-clause BSD style license.
