#!/usr/bin/env python
# encoding: utf-8
"""
createWorkableDB.py

Created by José Sánchez-Gallego on 6 Dec 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

import psycopg2
from psycopg2.extras import DictCursor
import os
import sys
import subprocess


USER = 'albireo'
PASSWORD = ''

SCHEMA = os.path.join(os.path.dirname(__file__), 'plateDB_Test.dump')
# SCHEMA = os.path.join(os.path.dirname(__file__), '../ignore/platedb.sql')
MANGADB_SCHEMA = os.path.join(os.path.dirname(__file__), 'mangaDB.sql')
OVERWRITE = True
DBNAME = 'plateDB_Test'
# DBNAME = 'apo_platedb'
ADD_LABELS = True
ADD_TILES = True


def _createConnection(database, user, password):
    conn = psycopg2.connect(database=database, host='localhost',
                            user=user, password=password)
    cur = conn.cursor(cursor_factory=DictCursor)
    return conn, cur


def addLabels(dbName):

    print 'Adding labels'

    conn, cur = _createConnection(dbName, USER, PASSWORD)

    tiles_status = [[0, 'Automatic'],
                    [1, 'Do Not Observe'],
                    [2, 'Override Complete'],
                    [3, 'Override Incomplete']]
    for tiles_status_row in tiles_status:
        cur.execute('INSERT INTO platedb.tile_status VALUES (%s, %s)',
                    tiles_status_row)

    plate_completion_status = [[0, 'Automatic'],
                               [1, 'Do Not Observe'],
                               [2, 'Force Complete'],
                               [3, 'Force Incomplete']]
    for plate_completion_status_row in plate_completion_status:
        cur.execute('INSERT INTO platedb.plate_completion_status ' +
                    'VALUES (%s, %s)',
                    plate_completion_status_row)

    plate_location = [[0, 'Design'],
                      [20, 'UW'],
                      [21, 'Storage'],
                      [22, 'Merrelli\'s_Apartment'],
                      [23, 'APO'],
                      [24, 'In Transit To Storage'],
                      [25, 'Recycled'],
                      [26, 'Given Away']]
    for plate_location_row in plate_location:
        cur.execute('INSERT INTO platedb.plate_location VALUES (%s, %s)',
                    plate_location_row)

    conn.commit()
    conn.close()


# Checks if the database exists and removes it if OVERWRITE is TRUE
conn, cur = _createConnection('postgres', 'postgres', 'postgres')
cur.execute('select datname from pg_database')
dbList = zip(*cur.fetchall())[0]

if DBNAME in dbList:
    if OVERWRITE:
        print 'Removing database {0}.'.format(DBNAME)
        cc = subprocess.Popen('dropdb {0}'.format(DBNAME), shell=True)
        cc.communicate()
    else:
        print 'Database {0} exists.'.format(DBNAME)
        sys.exit(1)

# Creates the database from a template
print 'Creating database'
cc = subprocess.Popen('createdb -T template0 {0}'.format(DBNAME),
                      shell=True)
cc.communicate()

# Restores the schema
print 'Populating database'
cc = subprocess.Popen('psql {0} < {1}'.format(DBNAME, SCHEMA), shell=True,
                      stdout=subprocess.PIPE)
cc.communicate()

# Adds labels and generic fields to the database
if ADD_LABELS:
    addLabels(DBNAME)

# Creates mangaDB schema
print 'Creating mangaDB schema'

cc = subprocess.Popen('psql {0} < {1}'.format(DBNAME, MANGADB_SCHEMA),
                      shell=True,
                      stdout=subprocess.PIPE)
cc.communicate()


if ADD_TILES:
    from Totoro.helpers import addFromTilingCatalogue
    print 'Adding tiles'
    samplePath = os.environ['MANGASAMPLE']
    addFromTilingCatalogue(samplePath, isSample=True)
