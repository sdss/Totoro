#!/usr/bin/env python
# encoding: utf-8
"""

utils.py

Created by José Sánchez-Gallego on 15 Oct 2015.
Licensed under a 3-clause BSD license.

Revision history:
    15 Oct 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division, print_function

import os
import re
# import inflect
import warnings

from SDSSconnect.exceptions import SDSSconnectUserWarning


classConversions = {
    'Mangadb_Sn2Values': 'Mangadb_SN2Values',
    'Mangadb_ExposureToDataCube': 'Mangadb_ExposureToData_cube',
    'Platedb_PlPlugmapM': 'Platedb_PlPlugMapM',
    'Platedb_BossSn2Threshold': 'Platedb_BossSN2Threshold',
    'Platedb_ProfMeasurement': 'Platedb_ProfilometryMeasurement',
    'Platedb_ProfTolerances': 'Platedb_ProfilometryTolerances'
}


def scalarRelationship(base, local_cls, referred_cls, constraint):
    """Creates relationship name removing the schema prefix."""

    referred_name = referred_cls.__name__
    if (str(referred_cls.__table__.schema) == str(local_cls.__table__.schema)):
        referred_name = referred_name.replace(
            str(referred_cls.__table__.schema).capitalize() + '_', '')

    uncamelised = re.sub(r'[A-Z]', lambda m: '_%s' % m.group(0).lower(), referred_name)[1:]

    return uncamelised


def cameliseClassname(base, tablename, table):
    """Produces a 'camelised' class name."""

    camelised = (str(table.schema).capitalize() + '_' + str(
        tablename[0].upper() + re.sub(r'_([a-z])', lambda mm: mm.group(1).upper(), tablename[1:])))

    if camelised in classConversions:
        return classConversions[camelised]

    return camelised


# _pluralizer = inflect.engine()
#
#
# def pluraliseCollection(base, local_cls, referred_cls, constraint):
#     """Produce an 'uncamelised', 'pluralised' class name."""
#
#     referred_name = referred_cls.__name__
#
#     if str(referred_cls.__table__.schema) == str(local_cls.__table__.schema):
#         referred_name = referred_name.replace(
#             str(referred_cls.__table__.schema).capitalize() + '_', '')
#
#     uncamelised = re.sub(r'[A-Z]',
#                          lambda m: "_%s" % m.group(0).lower(),
#                          referred_name)[1:]
#     pluralised = _pluralizer.plural(uncamelised)
#     pluralised = str(pluralised).replace('__', '_')
#
#     return pluralised


def modelRepr(self):
    """Creates a custom representation label for Model Classes."""

    className = str(self.__class__.__name__)
    if '_' in className:
        schema, table = className.split('_')
        schema = schema.capitalize()
        if schema[-2:] == 'db':
            schema = schema[0:-2] + 'DB'
    else:
        schema = ''
        table = className

    fullClassName = '{0}.{1}'.format(schema, table) if schema != '' else table

    reprList = ['pk={0}'.format(self.pk)]

    reprFile = os.path.join(os.path.dirname(__file__), '../etc/repr.dat')
    for line in open(reprFile):

        if line.strip().startswith('#') or line.strip() == '':
            continue

        try:
            reprSchema, reprTable, reprLabel, reprRef = line.strip().split()
            if (schema.lower() == reprSchema.lower() and table.lower() == reprTable.lower()):
                reprList.append('{0}={1}'.format(reprLabel, eval('self.{0}'.format(reprRef))))
        except Exception:
            warnings.warn('failed adding custom representation to {0}'.format(fullClassName),
                          SDSSconnectUserWarning)

    return '<{0}: {1}>'.format(fullClassName, ', '.join(reprList))


def nullifyRelationship(base, direction, return_fn, attrname, local_cls, referred_cls, **kwargs):

    return None


class ModelWrapper(object):

    def __init__(self, Base, prefix):
        """Creates a wrapper around Base.class names with a certain prefix."""

        for key in list(Base.classes.keys()):
            if prefix in key:
                exec('self.{0} = Base.classes.{1}'.format(key[len(prefix):], key))
                exec('self.{0}.__repr__ = modelRepr'.format(key[len(prefix):]))
