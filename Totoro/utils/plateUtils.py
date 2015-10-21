#!/usr/bin/env python
# encoding: utf-8
"""
plateUtils.py

Created by José Sánchez-Gallego on 7 Mar 2014.
Licensed under a 3-clause BSD license.

A collection of convenience functions to interact with plateDB.

"""

__ALL__ = ['isPlateCompleted', 'field2plate',
           'design2plate', 'getCompletedFields']


from ..exceptions import TotoroError
from ..plateDB.dataModel import Plate, MangaDB_Field_To_Plate
from .. import session
from ..core.defaults import R_SN2, B_SN2


def isPlateCompleted(locationIDs=None, plateIDs=None, designIDs=None):
    """Determines if a plate has been completed.

    This function accepts a list of loactionIDs, plateIDs or designIDs
    and determines if they have been completed. Either ``locationIDs``
    or ``plateIDs`` or ``designIDs`` must be specified with an integer
    or a list of integers. If more than one of the parameters is defined,
    an error is raised. The returned value is a bool or tuple of bools
    of the same size as the input. True indicates that the plate has
    been completed, False that it has not been observed. None values
    indicate that the input location, plate or design does not exist
    in plateDB.

    Parameters
    ----------
    locationIDs : int or list of ints, optional
        The locationID(s) to be checked.
    plateIDs : int or list of ints, optional
        The plateID(s) to be checked.
    designIDs : int or list of ints, optional
        The designID(s) to be checked.

    Result
    ------
    isCompleted : bool, None or tuple
        If the design, plate or location is in plateDB and has been
        completed, returns True; False otherwise. None indicates that
        the input value cannot be found in plateDB. If the input value
        is a list, a tuple of the same size is returned.

    """

    if locationIDs is None and plateIDs is None and designIDs is None:
        raise TotoroError('no input.')

    elif (locationIDs is not None and plateIDs is not None) or \
            (locationIDs is not None and designIDs is not None) or \
            (plateIDs is not None and designIDs is not None):
        raise TotoroError('only one type of input is permitted.')

    elif locationIDs is not None:
        inputs = locationIDs
        inputType = 'location'

    elif plateIDs is not None:
        inputs = plateIDs
        inputType = 'plate'

    elif designIDs is not None:
        inputs = designIDs
        inputType = 'design'

    if hasattr(inputs, '__getitem__'):
        isList = True
    else:
        isList = False
        inputs = [inputs]

    if inputType == 'location':
        plates = [field2plate(ii) for ii in inputs]
    elif inputType == 'design':
        plates = [design2plate(ii) for ii in inputs]
    else:
        plates = [plateid2platepk(ii) for ii in inputs]

    completed = []
    for plate in plates:

        if not hasattr(ii, '__getitem__'):
            completed.append(_isPlateCompleted(plate))
        else:
            completed.append([_isPlateCompleted(pp) for pp in plate])

    for ii in range(len(completed)):
        if hasattr(completed[ii], '__getitem__'):
            if None in completed[ii]:
                completed[ii] = None
            elif False in completed[ii]:
                completed[ii] = False
            else:
                completed[ii] = True

    if isList:
        return completed
    else:
        return completed[0]


def _isPlateCompleted(platePK):
    """Helper function to determine if a plate is completed."""

    if platePK is None:
        return None

    plate = session.query(Plate).filter(Plate.plate_pk == platePK)
    if plate.count() == 0:
        return None

    dataCubeQuery = plate[0].mangaDB_dataCube

    if len(dataCubeQuery) == 0:
        return None

    dataCube = dataCubeQuery[0]

    if dataCube.r1_sn2 >= R_SN2 and dataCube.r2_sn2 >= R_SN2 and \
            dataCube.b1_sn2 >= B_SN2 and dataCube.b2_sn2 >= B_SN2:
        return True
    else:
        return False


def field2plate(field):
    """Determines the plateID(s) for a field."""

    query = session.query(MangaDB_Field_To_Plate).filter(
        MangaDB_Field_To_Plate.field_pk == field)
    if query.count() == 0:
        return None
    else:
        return [qq.plate_pk for qq in query]


def design2plate(design):
    """Determines the plateID for a locationID."""

    query = session.query(Plate.pk).filter(Plate.design_pk == design)
    if query.count() == 0:
        return None
    else:
        return query.all()


def plateid2platepk(plateid):

    query = session.query(Plate.pk).filter(Plate.plate_id == plateid)
    if query.count() == 0:
        return None
    else:
        return query.all()


def getCompletedFields():

    fields = session.query(MangaDB_Field_To_Plate)

    if fields.count() == 0:
        return []
    else:
        completed = []
        for field in fields:
            if _isPlateCompleted(field.plate_pk) is True:
                completed.append(fields.location_id)
        return completed
