
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

CREATE TABLE field (pk SERIAL PRIMARY KEY NOT NULL,
    name TEXT, center_ra REAL, center_dec REAL,
    location_id INTEGER, expected_no_visits INTEGER,
    shared BOOLEAN DEFAULT FALSE, field_type_pk INTEGER,
    field_length_pk INTEGER);

CREATE TABLE field_type (pk SERIAL PRIMARY KEY NOT NULL,
    label TEXT);

CREATE TABLE field_length (pk SERIAL PRIMARY KEY NOT NULL,
    label TEXT);

CREATE TABLE data_cube (pk SERIAL PRIMARY KEY NOT NULL,
    plate_pk INTEGER, r1_sn2 REAL, r2_sn2 REAL, b1_sn2 REAL,
    b2_sn2 REAL);

CREATE TABLE field_to_plate (pk SERIAL PRIMARY KEY NOT NULL,
    field_pk INTEGER, plate_pk INTEGER);

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

ALTER TABLE ONLY field
    ADD CONSTRAINT field_type_fk FOREIGN KEY (field_type_pk)
    REFERENCES field_type(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY field
    ADD CONSTRAINT field_length_fk FOREIGN KEY (field_length_pk)
    REFERENCES field_length(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY data_cube
    ADD CONSTRAINT platedb_plate_fk FOREIGN KEY (plate_pk)
    REFERENCES platedb.plate(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY spectrum
    ADD CONSTRAINT data_cube_fk
    FOREIGN KEY (data_cube_pk) REFERENCES data_cube(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY field_to_plate
    ADD CONSTRAINT field_fk
    FOREIGN KEY (field_pk) REFERENCES field(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

ALTER TABLE ONLY field_to_plate
    ADD CONSTRAINT field_to_plate_plate_fk
    FOREIGN KEY (plate_pk) REFERENCES platedb.plate(pk)
    ON UPDATE CASCADE ON DELETE CASCADE;

INSERT INTO set_status VALUES (0, 'Incomplete'), (1, 'Excellent'), (2, 'Good'), (3, 'Poor');
INSERT INTO exposure_status VALUES (0, 'Good'), (1, 'Bad');
