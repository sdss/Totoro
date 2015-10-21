#!/usr/bin/env python
# encoding: utf-8
"""
defaults.py

Created by José Sánchez-Gallego on 7 Mar 2014.
Licensed under a 3-clause BSD license.

Default values for commonly used parameters.

"""


from .configuration import ConfigObject
import os

EXPTIME = ConfigObject('expTime', 15, 'Exposure time in minutes.')
EFFICIENCY = ConfigObject('efficiency', 0.75, 'Observation efficiency')
NDITHERS = ConfigObject('nDithers', 3, 'Number of dithers')

APO_LONGITUDE = 254.179722
LONGITUDE = ConfigObject('longitude', APO_LONGITUDE,
                         'The longitude of the observatory',
                         section='Observatory')

APO_LATITUDE = 32. + 46. / 60.
LATITUDE = ConfigObject('latitude', APO_LATITUDE,
                        'The latitude of the observatory.',
                        section='Observatory')

SURVEY = ConfigObject('survey', 'MaNGA', 'The survey name.')
SAMPLE = ConfigObject('sample_path', None, 'The path to the sample catalogue')

AVG_SEEING = ConfigObject(
    'avgSeeing', 1.0, 'The average seeing to use in the simulations.',
    section='Seeing')

MIN_SEEING_EXPOSURE = ConfigObject(
    'minSeeingExposure', 2.5,
    'Minimum average seeign for a exposure to be valid.',
    section='Seeing')

SEEING_POOR = ConfigObject(
    'seeingPoor', 2., 'Threshold above which the seeing is considered poor.')

SEEING_EXCELLENT = ConfigObject(
    'seeingExcellent', 1.5, 'Seeing below which the seeing is excellent.')

MAX_DIFF_SEEING_SET = ConfigObject(
    'maxDiffSeeingSet', 0.8,
    'Maximum difference between seeing in one set to be considered valid.',
    section='Seeing')

AVG_SN_RED = ConfigObject(
    'r_sn2_average', 7.5,
    'Average SN2 in the red band in a 15-minute exposure',
    section='SN')

AVG_SN_BLUE = ConfigObject(
    'b_sn2_average', 3.6,
    'Average SN2 in the red band in a 15-minute exposure',
    section='SN')

ALPHA_RED = ConfigObject(
    'alphaRed', 1.25, 'Airmass correction exponent, red band.')

ALPHA_BLUE = ConfigObject(
    'alphaBlue', 1., 'Airmass correction exponent, blue band.')

R_SN2 = ConfigObject('r_sn2', 60,
                     'Threshold SN2 for the R1 and R2 spectrographs.',
                     section='SN')

B_SN2 = ConfigObject('b_sn2', 27,
                     'Threshold SN2 for the B1 and B2 spectrographs.',
                     section='SN')

R_SN2_EXPOSURE = ConfigObject(
    'r_sn2_exposure', 3., 'Threshold SN2 for the R spectrograph per exposure',
    section='SN')

B_SN2_EXPOSURE = ConfigObject(
    'b_sn2_exposure', 3., 'Threshold SN2 for the B spectrograph per exposure',
    section='SN')

SN2_FACTOR_SET = ConfigObject(
    'sn2FactorSet', 2., section='SN')

DITHER_POSITIONS_NEEDED = ConfigObject(
    'ditherPositionsNeeded', ['N', 'S', 'E'],
    'Needed dither positions for a dither to be compete.')
NDITHERS = ConfigObject('nDithers', 3, 'Minimum number of exposures in a set.')
DITHER_POSITIONS = ConfigObject(
    'ditherPosition', ['C', 'E', 'S', 'N'], 'Accepted dither positions.')
EXPTYPES = ConfigObject('obsTypes', ['sci', 'cal'], 'Accepted exposure types.')


DEFAULT_PLAN_FILE = ConfigObject('defaultPlanFile',
                                 os.path.join(
                                     os.path.dirname(__file__),
                                     '../data/nightly.D.txt'),
                                 'The file with the observing plan',
                                 section='Scheduling')

OPTIMISED_PLAN_PATTERN = ConfigObject('optimisedPlanPattern',
                                      '{0}.jd',
                                      'The format for the optimised plan',
                                      section='Scheduling')

MJD_COLNAMES = ['MJD', 'APOGEE_0', 'APOGEE_1', 'MaNGA_0', 'MaNGA_1',
                'eBOSS_0', 'eBOSS_1']
AUTOSCHEDULER_COLNAMES = ['col{0}'.format(ii) for ii in range(1, 14)]
AUTOSCHEDULER_VALID_COLNAMES = ['col1', 'col4', 'col5', 'col8',
                                'col9', 'col12', 'col13']

INCREASE_MAPS = ConfigObject(
    'increaseMaps',
    os.path.join(
        os.path.dirname(__file__),
        '../data/IGincrease.fits'),
    'The g- and i-band increase maps.')
INCREASE_MAPS_FORMAT = ConfigObject('increaseMapsFormat', 'grid')

SKY_PRIORITY = ConfigObject('skyPriority', 4, section='Observing Logic')

COMPLETION_PRIORITY = ConfigObject('completionPriority', 6,
                                   section='Observing Logic')
