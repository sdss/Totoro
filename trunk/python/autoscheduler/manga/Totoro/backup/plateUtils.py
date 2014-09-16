#!/usr/bin/env python
# encoding: utf-8
"""
plateUtils.py

Created by José Sánchez-Gallego on 7 Mar 2014.
Licensed under a 3-clause BSD license.

A collection of convenience functions to interact with plateDB.

"""


# from ..plateDB.dataModel import (Plate, MangaDB_Field_To_Plate,
#                                  Plugging, MangaDB_Field)
# from ..plateDB import session
from ..core.defaults import LATITUDE
import numpy as np
from numbers import Real


# def field2plate(field_pk):
#     """Determines the plateID(s) for a field."""

#     query = session.query(MangaDB_Field_To_Plate).filter(
#         MangaDB_Field_To_Plate.field_pk == field_pk)
#     if query.count() == 0:
#         return None
#     else:
#         return [qq.plate_pk for qq in query]


# def design2plate(design):
#     """Determines the plateID for a locationID."""

#     query = session.query(Plate.pk).filter(Plate.design_pk == design)
#     if query.count() == 0:
#         return None
#     else:
#         return query.all()


# def plateid2platepk(plateid):

#     query = session.query(Plate.pk).filter(Plate.plate_id == plateid)
#     if query.count() == 0:
#         return None
#     else:
#         return query.all()


# def areFieldsDone(field_pks):
#     """Checks if a field is flagged good in plateDB.

#     This function checks all the plates drilled for a field and
#     all the pluggings for each plate and returns True if all
#     have been flagged as good.

#     """

#     isList = True
#     if not isinstance(field_pks, (list, tuple, np.ndarray)):
#         field_pks = [field_pks]
#         isList = False

#     qq = session.query(
#         MangaDB_Field.pk, Plugging.plugging_status_pk).join(
#         MangaDB_Field_To_Plate).join(Plate).join(Plugging).all()

#     qq = np.array(qq)

#     doneStatus = []
#     for ff in field_pks:

#         if len(qq) == 0:
#             doneStatus.append(False)
#             continue

#         plugging_statuses = qq[qq[:, 0] == ff][:, 1]

#         if 1 in plugging_statuses or 3 in plugging_statuses:
#             doneStatus.append(True)
#         else:
#             doneStatus.append(False)

#     if isList:
#         return doneStatus
#     else:
#         return doneStatus[0]


