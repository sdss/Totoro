#!/usr/bin/env python
# encoding: utf-8
"""
getMaNGAstats.py

Created by José Sánchez-Gallego on 27 Jan 2015.
Licensed under a 3-clause BSD license.

Revision history:
    27 Jan 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.internal.manga.Totoro.dbclasses.plate import getAll
from sdss.internal.manga.Totoro.logic.mangaLogic import checkExposure
from astropy import table
import numpy as np
import argparse
import sys
import os


def createWikiOutput(expTable, statTable, wiki):
    """Outputs the exposure table in Trac wiki format."""

    orphanExpPer = (np.sum(expTable['orphan_exposures']) /
                    np.sum(expTable['valid_exposures']) * 100.)
    invalidExpPer = (np.sum(expTable['invalid_exposures']) /
                     np.sum(expTable['exposures']) * 100.)

    for col in expTable.colnames:
        expTable.rename_column(
            col, '\'\'\'' + col.replace('_', ' ').capitalize() + '\'\'\'')
        statTable.rename_column(
            col, '\'\'\'' + col.replace('_', ' ').capitalize() + '\'\'\'')

    expTable.write(wiki, format='ascii.fixed_width', delimiter='  ||  ')

    unit = open(wiki, 'a')
    statTable.write(unit, format='ascii.fixed_width', delimiter='  ||  ')

    unit.write('\n\'\'\'Orphan exposures:\'\'\' {0:.1f}% [[BR]]'
               .format(orphanExpPer))
    unit.write('\n\'\'\'Invalid exposures:\'\'\' {0:.1f}%'
               .format(invalidExpPer))
    unit.write('\n')

    unit.close()

    # Reopens the file to remove the header from statTable
    lines = open(wiki, 'r').read().splitlines()

    unit = open(wiki, 'w')

    for nn, line in enumerate(lines):
        if nn == len(expTable) + 1:
            continue
        unit.write(line + '\n')

    unit.close()

    return


def getMaNGAstats(output=None, wiki=False, quiet=False):
    """Returns a table with some MaNGA statistics.

    This function returns a table (that can also be output to a file or in
    Trac Wiki format) with statistics of the exposures for complete plates.

    Keyword arguments
    -----------------
    output : string or None
        If output is a string, the table is saved to a file with that name in
        fixed width column format.

    wiki : False, None or string
        If False, not output in Trac Wiki format is produced. If a string, the
        wiki output is saved to that file. If None, the output is saved to a
        file with default filename "exposures_trac.dat"

    quiet : bool
        If True, the table is not output to stdout.

    Returns
    -------
    result : `astropy.table.Table`
        An `astropy.table.Table` in which each row is a plate and the columns
        are the different statistics being calculated. The table is also
        printed unless `quiet=True`.

    """

    allPlates = getAll(updateSets=False, silent=True, fullCheck=False)

    complete = [plate for plate in allPlates
                if plate.isComplete and plate.plate_id > 7400]

    names = ['plate_id', 'exposures', 'valid_exposures',
             'exposures_in_valid_sets', 'orphan_exposures',
             'invalid_exposures', 'invalid_SN2', 'invalid_seeing',
             'invalid_HA', 'invalid_other']
    expTable = table.Table(None, names=names, dtype=10 * [int])
    statTable = table.Table(None, names=names, dtype=['S10'] + 9 * [float])

    for plate in complete:

        plateID = str(plate.plate_id)
        exposures = plate.getMangadbExposures()
        nExposures = len(exposures)

        validExposures = []
        invalidExposures = []
        for exp in exposures:
            if checkExposure(exp.pk, parent='mangaDB', silent=True)[0] is True:
                validExposures.append(exp)
            else:
                invalidExposures.append(exp)

        nValidExposures = len(validExposures)
        nInvalidExposures = len(invalidExposures)

        exposuresInValidSets = 3 * len([ss for ss in plate.sets
                                        if ss.getQuality()[0] not in
                                        ['Incomplete', 'Bad', 'Unplugged']])

        orphanExposures = nValidExposures - exposuresInValidSets

        invalidSN2 = invalidSeeing = invalidHA = invalidOther = 0

        for exp in invalidExposures:
            pk = exp.pk
            status, errorCode = checkExposure(pk, parent='mangaDB',
                                              flag=False, force=True,
                                              silent=True)
            if errorCode == 3:
                invalidSeeing += 1
            elif errorCode == 4 or errorCode == 7:
                # errorCode=7 is low transparency but here we combine that
                # with low SN2.
                invalidSN2 += 1
            elif errorCode == 5:
                invalidHA += 1
            elif errorCode != 0 and errorCode != 10:
                invalidOther += 1

        expTable.add_row((plateID, nExposures, nValidExposures,
                          exposuresInValidSets, orphanExposures,
                          nInvalidExposures, invalidSN2, invalidSeeing,
                          invalidHA, invalidOther))

    statTable.add_row(['Total'] + [np.round(np.sum(expTable[col]), 0)
                                   for col in expTable.colnames[1:]])
    statTable.add_row(['Mean'] + [np.round(np.mean(expTable[col]), 1)
                                  for col in expTable.colnames[1:]])
    statTable.add_row(['Median'] + [np.round(np.median(expTable[col]), 0)
                                    for col in expTable.colnames[1:]])

    orphanExpPer = (np.sum(expTable['orphan_exposures']) /
                    np.sum(expTable['valid_exposures']) * 100.)
    invalidExpPer = (np.sum(expTable['invalid_exposures']) /
                     np.sum(expTable['exposures']) * 100.)

    statTableFmt = statTable.pformat()

    if not quiet:
        expTable.pprint()
        for row in statTableFmt[2:]:
            print(row)
        print('\nOrphan exposures: {0:.1f}%'.format(orphanExpPer))
        print('Invalid exposures: {0:.1f}%'.format(invalidExpPer))

    if output is not None:
        assert isinstance(output, basestring)
        expTable.write(output, format='ascii.fixed_width')

    if wiki is not False:
        if wiki is None:
            wiki = 'exposures_trac.dat'
        assert isinstance(wiki, basestring)
        createWikiOutput(expTable, statTable, wiki)

    return expTable


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                                     prog=os.path.basename(sys.argv[0]))
    parser.add_argument('-o', '--output', metavar='output', type=str,
                        help='The output file.')
    parser.add_argument('-w', '--wiki', dest='wiki', metavar='wiki', type=str,
                        default=False,
                        help='The file to which the Trac wiki data ' +
                        'will be saved.')
    parser.add_argument('-q', '--quiet', dest='quiet', action='store_true',
                        help='Disables stdout output.')

    args = parser.parse_args()

    getMaNGAstats(output=args.output, wiki=args.wiki, quiet=args.quiet)
