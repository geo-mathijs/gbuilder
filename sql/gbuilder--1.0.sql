CREATE OR REPLACE FUNCTION ConcaveHull(geometry[], alpha float8 DEFAULT 1.0)
    RETURNS geometry
    AS 'MODULE_PATHNAME', 'concave_hull'
    LANGUAGE C STABLE STRICT;

CREATE OR REPLACE FUNCTION concave_hull_finalfn(internal)
    RETURNS geometry
    AS 'MODULE_PATHNAME'
    LANGUAGE C;

CREATE AGGREGATE ConcaveHull(
    BASETYPE = geometry,
    SFUNC = pgis_geometry_accum_transfn,
    STYPE = internal,
    FINALFUNC = concave_hull_finalfn
);
