/*----------------------------------------------------------------
 *
 * src/gbuilder.h
 *
 * All custom data structs that are used in gbuilder
 *
 *----------------------------------------------------------------
 */

#ifndef _GBUILDER_H
#define _GBUILDER_H

#include <utils/array.h>
#include <utils/geo_decls.h>

/*
 * Postgis struct used for array creation. The struct is private
 * so we just copy it over as we only need it once.
 */
typedef struct {
    ArrayBuildState *a;
    Datum data;
} pgis_abs;

/*
 * Edge contains indices to a POINT2D and Triangle contains
 * indices to a POINT2D plus a circumcirlce. Polygons should be
 * able to contain inner rings and have no need for a bounding
 * box, so the default POLYGON type is not used here.
 */
typedef struct {
    int p1;
    int p2;
} Edge;

typedef struct {
    float8 x;
    float8 y;
    float8 r;
    int p1;
    int p2;
    int p3;
    bool flag;
} Triangle;

typedef struct {
    float8 x[2];
    float8 y[2];
    bool flag;
} hull_t;

typedef struct {
    Point *pt;
    int size;
    int rings;
    int *index;
} polygon_t;

/*
 * Headers for the different geometry structs
 */
typedef struct {
    Triangle *tri;
    int len;
    int size;
} tri_t;

typedef struct {
    polygon_t *poly;
    int len;
    int size;
} multipolygon_t;

#endif /* _GBUILDER_H */
