#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2019-10-07
# @Filename: ic342_lst_range.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import sys

from Totoro.db import getConnection


def load_lst_ranges(lst0, lst1):

    lst0 = float(lst0)
    lst1 = float(lst1)

    ic342_ra = 55.910158

    ha0 = lst0 * 15. - ic342_ra
    ha1 = lst1 * 15. - ic342_ra

    db = getConnection()
    session = db.Session()

    ic342_plates = session.query(db.mangaDB.Plate).join(
        db.plateDB.Plate, db.plateDB.PlateToSurvey, db.plateDB.Survey,
        db.plateDB.SurveyMode).filter(db.plateDB.Survey.label == 'MaNGA',
                                      db.plateDB.SurveyMode.label == 'MaNGA 10min').all()

    with session.begin():
        for plate in ic342_plates:
            plate.ha_min = ha0
            plate.ha_max = ha1
            plate.field_name = 'IC342'


if __name__ == '__main__':

    lst0, lst1 = sys.argv[1:3]
    load_lst_ranges(ha0, lst1)
