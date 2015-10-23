#!/usr/bin/env python
# encoding: utf-8
'''

mangaDB.py

Created by José Sánchez-Gallego on 10 Oct 2015.
Licensed under a 3-clause BSD license.

Revision history:
    10 Oct 2015 J. Sánchez-Gallego
      Initial version

'''

from __future__ import division
from __future__ import print_function
from sqlalchemy import MetaData
from sqlalchemy.orm import relationship, backref
from SDSSconnect.models.utils import cameliseClassname, generateRelationship


def construct_mangaDB(engine, Base):

    metadataMangaDB = MetaData(schema='mangadb')
    Base.metadata = metadataMangaDB

    Base.prepare(engine, reflect=True,
                 classname_for_table=cameliseClassname,
                 generate_relationship=generateRelationship)

    # Relationships
    Base.classes.Mangadb_Exposure.set = relationship(
        Base.classes.Mangadb_Set, backref='exposures')
    Base.classes.Mangadb_Exposure.status = relationship(
        Base.classes.Mangadb_ExposureStatus, backref='exposures')
    Base.classes.Mangadb_Exposure.platedbExposure = relationship(
        Base.classes.Platedb_Exposure, backref='mangadbExposure')
    Base.classes.Mangadb_Exposure.spectra = relationship(
        Base.classes.Mangadb_Spectrum, backref='exposures')
    Base.classes.Mangadb_Exposure.datacubes = relationship(
        Base.classes.Mangadb_DataCube, backref='exposures')

    Base.classes.Mangadb_ExposureToData_cube.exposure = relationship(
        Base.classes.Mangadb_Exposure, backref='exposuresToDatacubes')
    Base.classes.Mangadb_ExposureToData_cube.datacube = relationship(
        Base.classes.Mangadb_DataCube, backref='exposuresToDatacubes')

    Base.classes.Mangadb_DataCube.plate = relationship(
        Base.classes.Platedb_Plate, backref='dataCube')

    Base.classes.Mangadb_Spectrum.datacube = relationship(
        Base.classes.Mangadb_DataCube, backref='spectrum')

    Base.classes.Mangadb_Set.status = relationship(
        Base.classes.Mangadb_SetStatus, backref='sets')

    Base.classes.Mangadb_SN2Values.exposure = relationship(
        Base.classes.Mangadb_Exposure, backref='sn2values')

    Base.classes.Mangadb_Plate.platedbPlate = relationship(
        Base.classes.Platedb_Plate, backref=backref('mangadbPlate',
                                                    uselist=False))

    return
