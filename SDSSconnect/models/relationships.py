#!/usr/bin/env python
# encoding: utf-8
"""

relationships.py

Created by José Sánchez-Gallego on 15 Oct 2015.
Licensed under a 3-clause BSD license.

Revision history:
    15 Oct 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sqlalchemy.orm import relationship, backref


def createRelationships(Base):
    """Creates relationships."""

    # PlateDB
    Base.classes.Platedb_PlateRun.plates = relationship(
        Base.classes.Platedb_Plate,
        order_by=Base.classes.Platedb_Plate.plate_id, backref='platerun')

    Base.classes.Platedb_Plate.design = relationship(
        Base.classes.Platedb_Design,
        primaryjoin=(Base.classes.Platedb_Plate.design_pk ==
                     Base.classes.Platedb_Design.pk), backref='plate')
    Base.classes.Platedb_Plate.location = relationship(
        Base.classes.Platedb_PlateLocation, backref='plates')
    Base.classes.Platedb_Plate.surveys = relationship(
        Base.classes.Platedb_Survey,
        secondary=Base.classes.Platedb_PlateToSurvey.__table__,
        backref='plates')
    Base.classes.Platedb_Plate.pluggings = relationship(
        Base.classes.Platedb_Plugging,
        order_by=Base.classes.Platedb_Plugging.fscan_mjd, backref='plate')
    Base.classes.Platedb_Plate.currentSurveyMode = relationship(
        Base.classes.Platedb_SurveyMode,
        primaryjoin=(Base.classes.Platedb_Plate.current_survey_mode_pk ==
                     Base.classes.Platedb_SurveyMode.pk), backref='plates')
    Base.classes.Platedb_Plate.statuses = relationship(
        Base.classes.Platedb_PlateStatus,
        secondary=Base.classes.Platedb_PlateToPlateStatus.__table__,
        backref='plates')
    Base.classes.Platedb_Plate.completionStatus = relationship(
        Base.classes.Platedb_PlateCompletionStatus, backref='plates')
    Base.classes.Platedb_Plate.cmmMeasurements = relationship(
        Base.classes.Platedb_CmmMeas, backref='plate')

    Base.classes.Platedb_PlatePointing.plate = relationship(
        Base.classes.Platedb_Plate, backref='plate_pointings')
    Base.classes.Platedb_PlatePointing.pointing = relationship(
        Base.classes.Platedb_Pointing, backref='plate_pointings')
    Base.classes.Platedb_PlatePointing.observations = relationship(
        Base.classes.Platedb_Observation,
        order_by=Base.classes.Platedb_Observation.mjd.desc(),
        backref='plate_pointing')

    Base.classes.Platedb_PlateCompletionStatusHistory.plate = relationship(
        Base.classes.Platedb_Plate, backref='completionStatusHistory')
    Base.classes.Platedb_PlateCompletionStatusHistory.completionStatus = \
        relationship(Base.classes.Platedb_PlateCompletionStatus,
                     backref='completionStatusHistory')

    Base.classes.Platedb_Tile.plates = relationship(
        Base.classes.Platedb_Plate,
        order_by=Base.classes.Platedb_Plate.plate_id, backref='tile')
    Base.classes.Platedb_Tile.status = relationship(
        Base.classes.Platedb_TileStatus, backref='tiles')

    Base.classes.Platedb_TileStatusHistory.tile = relationship(
        Base.classes.Platedb_Tile, backref='statusHistory')
    Base.classes.Platedb_TileStatusHistory.status = relationship(
        Base.classes.Platedb_TileStatus, backref='statusHistory')

    Base.classes.Platedb_Design.pointings = relationship(
        Base.classes.Platedb_Pointing, backref='design')
    Base.classes.Platedb_Design.values = relationship(
        Base.classes.Platedb_DesignValue, backref='design')
    Base.classes.Platedb_Design.inputs = relationship(
        Base.classes.Platedb_PlateInput, backref='design')

    Base.classes.Platedb_DesignValue.field = relationship(
        Base.classes.Platedb_DesignField, backref='design_values')

    Base.classes.Platedb_Plugging.cartridge = relationship(
        Base.classes.Platedb_Cartridge, backref='pluggings')
    Base.classes.Platedb_Plugging.plplugmapm = relationship(
        Base.classes.Platedb_PlPlugMapM, backref='plugging')
    Base.classes.Platedb_Plugging.instruments = relationship(
        Base.classes.Platedb_Instrument,
        secondary=Base.classes.Platedb_PluggingToInstrument.__table__,
        backref='pluggings')
    Base.classes.Platedb_Plugging.observations = relationship(
        Base.classes.Platedb_Observation, backref='plugging')
    Base.classes.Platedb_Plugging.activePlugging = relationship(
        Base.classes.Platedb_ActivePlugging, backref='plugging')
    Base.classes.Platedb_Plugging.status = relationship(
        Base.classes.Platedb_PluggingStatus, backref='pluggings')

    Base.classes.Platedb_Observation.status = relationship(
        Base.classes.Platedb_ObservationStatus, backref='observations')
    Base.classes.Platedb_Observation.exposures = relationship(
        Base.classes.Platedb_Exposure, backref='observation',
        order_by=(Base.classes.Platedb_Exposure.start_time,
                  Base.classes.Platedb_Exposure.exposure_no))

    Base.classes.Platedb_Exposure.camera = relationship(
        Base.classes.Platedb_Camera, backref='exposures')
    Base.classes.Platedb_Exposure.survey = relationship(
        Base.classes.Platedb_Survey, backref='exposures')
    Base.classes.Platedb_Exposure.flavor = relationship(
        Base.classes.Platedb_ExposureFlavor, backref='exposures')
    Base.classes.Platedb_Exposure.status = relationship(
        Base.classes.Platedb_ExposureStatus, backref='exposures')
    Base.classes.Platedb_Exposure.headerValues = relationship(
        Base.classes.Platedb_ExposureHeaderValue,
        order_by=Base.classes.Platedb_ExposureHeaderValue.index,
        backref='exposure')

    Base.classes.Platedb_ExposureHeaderValue.header = relationship(
        Base.classes.Platedb_ExposureHeaderKeyword, backref='headerValues')

    Base.classes.Platedb_Exposure.surveyMode = relationship(
        Base.classes.Platedb_SurveyMode, backref='exposures')

    Base.classes.Platedb_Camera.instrument = relationship(
        Base.classes.Platedb_Instrument, backref='cameras')

    Base.classes.Platedb_CameraFrame.camera = relationship(
        Base.classes.Platedb_Camera, backref='cameraFrames')
    Base.classes.Platedb_CameraFrame.exposure = relationship(
        Base.classes.Platedb_Exposure, backref='cameraFrames')

    Base.classes.Platedb_Gprobe.cartridge = relationship(
        Base.classes.Platedb_Cartridge, backref='gprobes')

    Base.classes.Platedb_BossPluggingInfo.plugging = relationship(
        Base.classes.Platedb_Plugging, backref='bossPluggingInfo')

    Base.classes.Platedb_BossSN2Threshold.camera = relationship(
        Base.classes.Platedb_Camera, backref='bossSN2Threshold')

    Base.classes.Platedb_Profilometry.plugging = relationship(
        Base.classes.Platedb_Plugging, backref='profilometries')
    Base.classes.Platedb_Profilometry.measurements = relationship(
        Base.classes.Platedb_ProfilometryMeasurement, backref='profilometry',
        order_by=Base.classes.Platedb_ProfilometryMeasurement.number,
        cascade='all, delete, delete-orphan')
    Base.classes.Platedb_Profilometry.tolerances = relationship(
        Base.classes.Platedb_ProfilometryTolerances, backref='profilometry')

    Base.classes.Platedb_ProfilometryTolerances.survey = relationship(
        Base.classes.Platedb_Survey, backref='profilometry_tolerances')

    Base.classes.Platedb_PlateHolesFile.plate = relationship(
        Base.classes.Platedb_Plate, backref='plateHolesFile')

    Base.classes.Platedb_PlPlugMapM.fibers = relationship(
        Base.classes.Platedb_Fiber, backref='plPlugMapM')

    Base.classes.Platedb_Fiber.plateHoles = relationship(
        Base.classes.Platedb_PlateHole, backref='fiber')

    Base.classes.Platedb_PlateHole.plateHoleType = relationship(
        Base.classes.Platedb_PlateHoleType, backref='plateHole')
    Base.classes.Platedb_PlateHole.plateHolesFile = relationship(
        Base.classes.Platedb_PlateHolesFile, backref='plateHole')
    Base.classes.Platedb_PlateHole.objectType = relationship(
        Base.classes.Platedb_ObjectType, backref='plateHole')

    # MangaDB
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
