
CREATE SCHEMA mangadb;

SET search_path TO mangadb;

CREATE TABLE exposure (pk SERIAL PRIMARY KEY NOT NULL,
    r1_sn2 REAL, r2_sn2 REAL, b1_sn2 REAL, b2_sn2 REAL,
    seeing REAL, transparency REAL, dither_number INTEGER,
    dither_ra REAL, dither_dec REAL, ha REAL,
    dither_position CHAR[1], comment TEXT, set_pk INTEGER,
    exposure_status_pk INTEGER, data_cube_pk INTEGER,
    platedb_exposure_pk INTEGER);

CREATE TABLE set (pk SERIAL PRIMARY KEY NOT NULL,
    comment TEXT, name TEXT, set_status_pk INTEGER);

CREATE TABLE set_status (pk SERIAL PRIMARY KEY NOT NULL,
    label TEXT);

CREATE TABLE exposure_status (pk SERIAL PRIMARY KEY NOT NULL,
    label TEXT);

CREATE TABLE plate (pk SERIAL PRIMARY KEY NOT NULL,
    plate_status_pk INTEGER, platedb_plate_pk INTEGER,
    tile_pk INTEGER);

CREATE TABLE plate_status (pk SERIAL PRIMARY KEY NOT NULL, label TEXT);

CREATE TABLE tile (pk SERIAL PRIMARY KEY NOT NULL,
    id INTEGER, ra_centre REAL, dec_centre REAL,
    name TEXT, priority INTEGER DEFAULT 4, platedb_tile_pk INTEGER);

CREATE TABLE survey_mode (pk SERIAL PRIMARY KEY NOT NULL,
    label TEXT);

CREATE TABLE data_cube (pk SERIAL PRIMARY KEY NOT NULL,
    plate_pk INTEGER, r1_sn2 REAL, r2_sn2 REAL, b1_sn2 REAL,
    b2_sn2 REAL);

CREATE TABLE spectrum (pk SERIAL PRIMARY KEY NOT NULL,
    data_cube_pk INTEGER, fiber INTEGER, exposure REAL, ifu_no INTEGER);

ALTER TABLE ONLY exposure
    ADD CONSTRAINT set_fk FOREIGN KEY (set_pk) REFERENCES set(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY exposure
    ADD CONSTRAINT exposure_status_fk
    FOREIGN KEY (exposure_status_pk) REFERENCES exposure_status(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY exposure
    ADD CONSTRAINT platedb_exposure_fk
    FOREIGN KEY (platedb_exposure_pk) REFERENCES platedb.exposure(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY exposure
    ADD CONSTRAINT data_cube_fk
    FOREIGN KEY (data_cube_pk) REFERENCES data_cube(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY set
    ADD CONSTRAINT set_status_fk FOREIGN KEY (set_status_pk)
    REFERENCES set_status(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY plate
    ADD CONSTRAINT plate_status_fk FOREIGN KEY (plate_status_pk)
    REFERENCES plate_status(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY plate
    ADD CONSTRAINT tile_fk FOREIGN KEY (tile_pk)
    REFERENCES tile(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY plate
    ADD CONSTRAINT platedb_plate_fk FOREIGN KEY (platedb_plate_pk)
    REFERENCES platedb.plate(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY tile
    ADD CONSTRAINT platedb_tile_fk FOREIGN KEY (platedb_tile_pk)
    REFERENCES platedb.tile(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY data_cube
    ADD CONSTRAINT plate_fk FOREIGN KEY (plate_pk)
    REFERENCES plate(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY spectrum
    ADD CONSTRAINT data_cube_fk
    FOREIGN KEY (data_cube_pk) REFERENCES data_cube(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

INSERT INTO set_status VALUES (0, 'Incomplete'), (1, 'Excellent'), (2, 'Good'), (3, 'Poor');
INSERT INTO exposure_status VALUES (0, 'Good'), (1, 'Bad');
INSERT INTO plate_status VALUES (0, 'Incomplete'), (1, 'Complete'),
    (2, 'Force Incomplete'), (3, 'Force Complete');
INSERT INTO survey_mode VALUES (0, 'MaNGA Dither'), (1, 'APOGEE Lead'), (2, 'Testing');
