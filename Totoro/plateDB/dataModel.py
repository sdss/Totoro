#!/usr/bin/env python
# encoding: utf-8
"""
dataModel.py

Created by José Sánchez-Gallego on 31 Oct 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

Classes for sqlalchemy describing the PlateDB data model.

"""

import sys
from sqlalchemy.orm import relationship, backref

try:
    from .connection import DatabaseConnection
    db = DatabaseConnection()
except ImportError:
    print('Could not find DatabaseConnection\n'
          'Did you set it up before running?')
    sys.exit(1)

Base = db.Base
session = db.session


#####################  PLATEDB CLASSES #################

# The data model for plateDB is incomplete and only covers
# the needs of Totoro.

class Exposure(Base):
    __tablename__ = 'exposure'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    observation = relationship('Observation',
                               backref='exposure')
    survey = relationship('Survey', backref='exposures')
    exposureStatus = relationship('ExposureStatus',
                                  backref='exposures')
    exposureFlavor = relationship('ExposureFlavor',
                                  backref='exposures')

    def __repr__(self):
        return '<Exposure {0} (pk={1})>'.format(self.survey.label, self.pk)


class Survey(Base):
    __tablename__ = 'survey'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    def __repr__(self):
        return '<Survey {0} (pk={1})>'.format(self.label, self.pk)


class ExposureStatus(Base):
    __tablename__ = 'exposure_status'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    def __repr__(self):
        return '<ExposureStatus {0} (pk={1})>'.format(self.label, self.pk)


class ExposureFlavor(Base):
    __tablename__ = 'exposure_flavor'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    def __repr__(self):
        return '<ExposureFlavor {0} (pk={1})>'.format(self.label, self.pk)


class CameraFrame(Base):
    __tablename__ = 'camera_frame'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    camera = relationship('Camera', backref='cameraFrames')
    exposure = relationship('Exposure',
                            backref=backref('cameraFrames'))

    def __repr__(self):
        return '<CameraFrame Camera={0} (pk={1}, exposure_pk={2})>'.format(
            self.camera.label, self.pk, self.exposure.pk)


class Camera(Base):
    __tablename__ = 'camera'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    instrument = relationship('Instrument', backref='camera')

    def __repr__(self):
        return '<Camera {0} (pk={1}, instrument={2})>'.format(
            self.label, self.pk, self.instrument.label)


class Instrument(Base):
    __tablename__ = 'instrument'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    def __repr__(self):
        return '<Instrument {0} (pk={1})>'.format(self.label, self.pk)


class Observation(Base):
    __tablename__ = 'observation'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}
    platePointing = relationship('PlatePointing', backref='observation')

    def __repr__(self):
        return '<Observation (pk={0})>'.format(self.pk)


class PlatePointing(Base):
    __tablename__ = 'plate_pointing'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}
    pointing = relationship('Pointing', backref='platePointing')

    def __repr__(self):
        return '<PlatePointing (pk={0})>'.format(self.pk)


class Pointing(Base):
    __tablename__ = 'pointing'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    def __repr__(self):
        return '<Pointing (pk={0})>'.format(self.pk)


class Plate(Base):
    __tablename__ = 'plate'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    plate_run = relationship('PlateRun', backref='plates')
    tile = relationship('Tile', backref='plates')

    def __repr__(self):
        return '<Plate (pk={0})>'.format(
            self.pk)


class PlateRun(Base):
    __tablename__ = 'plate_run'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    def __repr__(self):
        return '<PlateRun (pk={0}, label={1}, year={2})>'.format(
            self.pk, self.label, self.year)


class Design(Base):
    __tablename__ = 'design'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    def __repr__(self):
        return '<Design (pk={0}, comment={1}>'.format(
            self.pk, self.comment)


class Tile(Base):
    __tablename__ = 'tile'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    def __repr__(self):
        return '<Tile (pk={0}, id={1})>'.format(
            self.pk, self.id)


class Plugging(Base):
    __tablename__ = 'plugging'
    __table_args__ = {'autoload': True, 'schema': 'platedb'}

    plate = relationship('Plate', backref='pluggings')

    def __repr__(self):
        return '<Plugging (pk={0}, plate_id={1})>'.format(
            self.pk, self.plate.plate_id)



###################### MANGADB CLASSES #####################


class MangaDB_Set(Base):
    __tablename__ = 'set'
    __table_args__ = {'autoload': True, 'schema': 'mangadb'}

    setStatus = relationship('MangaDB_SetStatus', backref='mangadbSets')

    def __repr__(self):
        return '<MangaDB Set (pk={0}, name={1})>'.format(
            self.pk, self.name)


class MangaDB_SetStatus(Base):
    __tablename__ = 'set_status'
    __table_args__ = {'autoload': True, 'schema': 'mangadb'}

    def __repr__(self):
        return '<MangaDB Set_Status (pk={0}, label={1})>'.format(
            self.pk, self.label)


class MangaDB_ExposureStatus(Base):
    __tablename__ = 'exposure_status'
    __table_args__ = {'autoload': True, 'schema': 'mangadb'}

    def __repr__(self):
        return '<MangaDB Exposure_Status (pk={0}, label={1})>'.format(
            self.pk, self.label)


class MangaDB_Spectrum(Base):
    __tablename__ = 'spectrum'
    __table_args__ = {'autoload': True, 'schema': 'mangadb'}

    dataCube = relationship('MangaDB_DataCube', backref='spectra')

    def __repr__(self):
        return '<MangaDB Spectrum (pk={0})>'.format(self.pk)


class MangaDB_DataCube(Base):
    __tablename__ = 'data_cube'
    __table_args__ = {'autoload': True, 'schema': 'mangadb'}

    plate = relationship('Plate', backref='mangaDB_dataCube')

    def __repr__(self):
        return '<MangaDB Data_Cube (pk={0})>'.format(self.pk)


class MangaDB_Exposure(Base):
    __tablename__ = 'exposure'
    __table_args__ = {'autoload': True, 'schema': 'mangadb'}

    set = relationship('MangaDB_Set', backref='mangadbExposures')
    exposureStatus = relationship('MangaDB_ExposureStatus',
                                  backref='mangadbExposures')
    dataCube = relationship('MangaDB_DataCube', backref='mangadbExposures')
    platedbExposure = relationship('Exposure', backref='mangadbExposure')

    def __repr__(self):
        return '<MangaDB Exposure (pk={0})>'.format(self.pk)


class MangaDB_Field(Base):
    __tablename__ = 'field'
    __table_args__ = {'autoload': True, 'schema': 'mangadb'}

    def __repr__(self):
        return '<MangaDB Field (pk={0})>'.format(self.pk)


class MangaDB_Field_To_Plate(Base):
    __tablename__ = 'field_to_plate'
    __table_args__ = {'autoload': True, 'schema': 'mangadb'}

    plate = relationship('Plate')
    field = relationship('MangaDB_Field')

    def __repr__(self):
        return '<MangaDB Field_To_Plate (pk={0})>'.format(self.pk)
