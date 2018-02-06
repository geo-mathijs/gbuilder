# gbuilder
Extension for PostgreSQL that adds various geometry builders. Only supports the creation of a concave hull at the moment, but more functionality will be added over time. The extension is compiled against liblwgeom, so that it can use the PostGIS geometry types directly.

## Dependencies
- PostgreSQL >= 9.6
- PostGIS >= 2.3
- [Fast Robust Predicates](https://www.cs.cmu.edu/~quake/robust.html)

## Installation
Build and install the extension:
```
make
sudo make install
```

And enable it in PostgreSQL:
```
CREATE EXTENSION postgis;
CREATE EXTENSION gbuilder;
```

## ConcaveHull
Construct a concave hull from a set of points. There is also a non-aggregate function available that accepts a custom alpha value.
```
SELECT ConcaveHull(geom)
FROM (
    SELECT (ST_DumpPoints(ST_GeneratePoints(ST_Buffer(ST_GeomFromText('LINESTRING(4.48187 51.94380,4.44632 51.90386,4.52564 51.91557)', 4326), 0.001, 'endcap=round join=round'), 100))).geom
) x;
```
