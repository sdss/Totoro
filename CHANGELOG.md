# Change Log

## [2.1.4] - unreleased

### Changed
- For MaNGa time, get only plates with `MaNGA dither` or `MaNGA 10min`. This prevents selecting plugged plates with `MaNGA Globular` survey mode.

### Fixed
- Fix check of incomplete S/N in a new exposure.
- Import `factorial` from `scipy.special`.


## [2.1.3] - 2019-04-18

### Added
- Two custom regions to the footprint.
- `PETUNIA_WARNING` log level that can be used to output warnings to Petunia.

### Fixed
- A bug in `getPlatesInFootprint` that would make `coords=True` return incorrect results.


## [2.1.2] - 2019-04-03

### Fixed
- Logger not imported in deprecated script ``rearrangeSets.py``.


## [2.1.1] - 2019-01-28

### Fixed
- Problem with logging important-level messages.


## [2.1.0] - 2018-12-13

### Changed
- Use 9mm for the MaNGA ferrule size.

### Refactored
- Uses new logger file and `log.warning` instead of `warnings.warn`.

### Fixed
- Does not remove auto-reflected relationships in `DatabaseConnection` since it breaks with SQLA>=1.3. That should not affect anything since they were deleted mostly to keep the namespace tidy.


## [2.0.0] - 2018-06-20

### Added
- All scripts are now consolidated under the `totoro` CLI. `mangaPluggingRequest`, `overrideSet`, and `rearrangeSets` are still available but will be deprecated in the future. Also, cleaned some old, deprecated scripts that had been in the `bin` directory for a while.

### Refactored
- Totoro is now Python 2/3 compatible after applying futurize to all files.
- Applied `isort`, `yapf`, and `unify` to all files.
- Removed dependencies from ``sdss_python_module`` by moving the necessary files to Totoro.
- Totoro now uses `pydl.pydlutils.yanny` for all Yanny files reading.

### Changed
- Changed colour of backup plates in log.
- In `Plugger`, a replug now only uses the cart of its previous plugging if the cart is empty or has a complete or not started MaNGA plate.
- The preferred order of available carts passed to the master autoscheduler for APOGEE time is now reversed. This will make APOGEE use the carts with more broken fibres first, leaving the best ones for MaNGA-led plates.
- Updated the footprint for ALFALFA NGC.

### Fixes
- In `footprint`, copies the list of vertices of the patch to plot before applying the origin value (e.g., for Mollweide projection). This prevents incorrect wrapping of regions if the function is called more than once.

[View changes](https://github.com/sdss/Totoro/compare/1.8.4...2.0.0)


## [1.8.4] - 2018-05-13

### Changed
- Renamed `Totoro` CLI to `totoro`.
- Improved logger formatting in `totoro` CLI.

### Fixed
- A problem with the logger not parsing attributes passed to the record.


## [1.8.3] - 2018-05-11

### Added
- Added a `plugging` command to `Totoro` that replicates the `mangaPluggingRequest` script.


## [1.8.2] - 2018-05-11

### Added
- Added warning when `--no-backup` is passed to `Totoro simulate`.
- Added `--version` option to `Totoro`.


## [1.8.1] - 2018-05-10

### Added
- Option `rejectBackup` in `Planner` to remove backup plates from the list of plates to schedule. Default is `False`.
- Option `excludeStarted` in `Planner` to remove plates that have been started (have valid exposures). Default is `False`.
- A command `Totoro` for CLI interaction with Totoro. At this time it only has a `simulate` subcommand.


## [1.8.0] - 2017-12-20

### Added
- Option `completionThreshold` in configuration to define the minimum completion needed to mark a plate as complete. Set to 0.985 by default
- `utils.isPlateComplete` now accepts a `write_apocomplete` option that, if the plate is complete, write the APOcomplete file to the corresponding directory in mangacore.
- `utils.isPlateComplete` now accepts a `mark_complete` option that, if the plate is complete, sets the plugging status to ``Complete``.

### Fixed
- Fixed a bug that would mark Unplugged sets as complete.
- Quick fix for the case when a plate cannot be allocated to a cart. This causes `Plugger._getCart` to return `None`. The current fix just ignores the plate, which can lead to unallocated time.
- A bug when defining the `field_name` in `checkExposure` that prevented exposure files from being restored.


## [1.7.2] - 2017-11-06

### Added
- A new section, `config['specialPrograms']` that allows to override SN2 completion thresholds for special plates. Implemented initially for IC342.
- Special program plates can be completed after a number of good sets.
- Added columns `ha_min`, `ha_max`, and `field_name` to `mangadb.Plate`. `Plate.getHA()` will use `ha_min, ha_max` from the database, if they are set. The `field_name` column allows to identify a plate as part of a special program.
- Added support for new survey mode `MaNGA 10min`.


## [1.7.1] - 2017-06-29

### Added
- `Plate.get_mastar()` method that returns MaStar plates (i.e., plates that are APOGEE2-MaNGA with APOGEE as the lead survey).
- Plugger avoids cart 2 for plates for which the holes are too close.

### Changed
- Reduced number of permutationLimitIncomplete sin it was causing timeout problems with Petunia.


## [1.7.0] - 2016-06-09

### Changed
- Improved efficiency in loading many plates by optimising queries and
lazy loading all information during the initialisation of Plate and Exposure.


## [1.6.2] - 2016-04-07

### Fixed
- A bug that gives higher priority to not started vs started plates when they
both are completed during the plugging/planning simulation


## [1.6.1] - ???

### Changed
- ???


## [1.6] - ???

### Changed
- ???


## [1.5.4] - ???

### Changed
- Updated default goodWeatherFraction to 0.5.


## [1.5.3] - 2016-01-19

### Changed
- Modified simulation efficiencies to 0.755 (to avoid problems with three exps.
not making the 1-hour window).


## [1.5.2] - 2016-01-18

### Added
- Initial buffer for plugger. If the observing block is in the second part of the night, the scheduling starts defaults.plugger.initialBufferMin before the official schedule time. Issues a series of warnings. If Plugger is called with the keyword useInitialBuffer=False, the initial buffer is not applied. Right now all plugger tests use this option to avoid changing the results. This should be fixed eventually.


## [1.5.1] - 2016-01-12

### Fixed
- Fixed a bad keyword in the set rearrangement script that made it fail.

### Refactored
- Improved get plates function.


## [1.5] - 2015-12-12

### Fixed
- Plate.addMockExposure only rearranges exposures in complete sets if the number of such exposures in the plate is <= 5. Otherwise, it can cause serious delays for plugging.

### Refactored
- Major refactoring and documenting of plate_utils, scheduler_utils, planner, and plugger.
- Most classes and functions in dbclasses can now be also accessed from the root of the module.


## [1.2.1] - 2015-11-13

### Added
- runTests.py now restores, if possible, the test DB before running. If the script is called with '--no-restore' it will skip the restoration.

### Fixed
- Fixed a bug in the plugging logic that would chose a plate for a single exposure when another started plate was available for that same LST window, but had lower increase in SN2. The code now looks for plates that already have signal (either real or simulated) and prioritises those.


## [1.2.0] - 2015-11-10

### Changed
- When adding new exposures found in the DB, Totoro won't try to rearrange the exposures in incomplete sets to optimise SN2 and patch dithers. The rearrangement still happens for simulated exposures (plugger and planner).


## [1.1.2] - 2015-11-09

### Changed
- Adds Totoro to the PYTHONPATH when running the rearrangeSets scripts. This fixes a problem when rearrangeSets is called from Petunia and Totoro is not in the PYTHONPATH.


## [1.1.1] - 2015-11-09

### Fixed
- Bug in set rearrangement when one or more of the exposures is invalid


## [1.1.0] - 2015-11-08

### Changed
- Seeing for simulated exposures changed from 1 to 1.5 arcsec and set as a configurable option in defaults.simulation.seeing.

### Fixed
- Plugger: fixes a bug when a plugged plate cannot be observed for an integer number of sets. In that case, the orphaned exposures are removed and the remaining time is to be observed with a different plate. However, if there is no other plate that can be used to observe at least one whole set, the original, already plugged plate must be favoured.


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
- Plugger: offline carts are now given higher cart_order. This is because we don't want APOGEE to use MaNGA offline carts for co-designed plates, if possible. When the plugger is running for a MaNGA night offline carts are given the highest priority after all scheduled carts. For non-MaNGA nights, offline carts are given higher priority than any other cart except carts with plates with incomplete sets.
- Module file

### Changed
- Totoro moved to its own repo.
- scheduler_utils: plates are sorted by completion so that, if several plates complete with the same number of exposures, the one with highest completion gets chosen.
- Now uses SDSSconnect for connecting to the DB. Requires passwords to be defined in .pgpass.
- DB classes now don't inherit from model classes. Instead, the DB query object is recorded as object.\_dbObject. DB attributes can still be accessed from the Totoro object via a custom `__getattr__` method.
