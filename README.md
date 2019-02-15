# gbuilder
Extension for PostgreSQL that adds various geometry builders. Only supports the creation of a concave hull at the moment. The extension is compiled against liblwgeom, so that it can use the PostGIS geometry types directly.

## Dependencies
- PostgreSQL >= 9.6
- PostGIS >= 2.5
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
Construct a concave hull from a set of points. There is also a non-aggregate function available that accepts a custom alpha value. It is useful to note that this concave hull algorithm is vastly superior to the internal PostGIS version in both speed and precision.
```
SELECT ConcaveHull(geom)
FROM planet_osm_point p
WHERE p."public_transport" = 'stop_position'
AND p.geom && ST_Buffer(ST_GeomFromText('POINT(4.48187 51.94380)', 4326), 0.05);
```
