#include <postgres.h>
#include <catalog/pg_type.h>
#include <fmgr.h>
#include <funcapi.h>
#include <liblwgeom.h>
#include <utils/array.h>
#include "gbuilder.h"

PG_MODULE_MAGIC;

Datum concave_hull(PG_FUNCTION_ARGS);
Datum concave_hull_finalfn(PG_FUNCTION_ARGS);

PG_FUNCTION_INFO_V1(concave_hull);
Datum concave_hull(PG_FUNCTION_ARGS) {
    ArrayIterator iterator;
    ArrayType *array;
    Datum value;
    GSERIALIZED *result;
    LWMPOLY * lwmpoly;
    LWPOINT *lwpoint;
    LWPOLY *lwpoly;
    POINT2D *points;
    POINT4D p;
    bool isnull;
    float8 alpha = 1.;
    int i, j, k;
    int nelems;
    int npoints;
    int srid = 4326;
    multipolygon_t *multi;
    size_t ret_size = 0;
    tri_t *tri;

    array = PG_GETARG_ARRAYTYPE_P(0);
    nelems = ArrayGetNItems(ARR_NDIM(array), ARR_DIMS(array));

    if (nelems == 0)
        PG_RETURN_NULL();

    /* Alloc 3 extra in advance */
    points = palloc((sizeof(POINT2D) * nelems) + 3);
    npoints = 0;

    iterator = array_create_iterator(array, 0, NULL);
    while (array_iterate(iterator, &value, &isnull)) {
        GSERIALIZED *pglwgeom;

        if (isnull)
            continue;

        pglwgeom = (GSERIALIZED *) DatumGetPointer(value);
        if (gserialized_get_type(pglwgeom) != POINTTYPE)
            continue;

        lwpoint = lwgeom_as_lwpoint(lwgeom_from_gserialized(pglwgeom));
        memcpy(&(points[npoints++]), getPoint_internal(lwpoint->point, 0), ptarray_point_size(lwpoint->point));
        lwgeom_free(lwpoint_as_lwgeom(lwpoint));
    }
    array_free_iterator(iterator);

    if (npoints == 0)
        PG_RETURN_NULL();

    multi = (multipolygon_t *) palloc0(sizeof(multipolygon_t));
    tri = (tri_t *) palloc0(sizeof(tri_t));

    /* Call the robust predicate library to initialize the epsilon */
    exactinit();

    triangulate(points, npoints, tri, srid);

    if (tri->len > 4)
        construct_concave_hull(points, npoints, tri, multi, alpha);

    lwmpoly = lwmpoly_construct_empty(srid, 0, 0);

    for (i = 0; i < multi->len; i++) {
        POINTARRAY **pas = palloc(sizeof(POINTARRAY *) * multi->poly[i].rings);

        for (j = 0; j < multi->poly[i].rings; j++) {
            pas[j] = ptarray_construct_empty(0, 0, 0);

            for (k = multi->poly[i].index[j]; k < multi->poly[i].index[j+1]; k++) {
                p.x = multi->poly[i].pt[k].x;
                p.y = multi->poly[i].pt[k].y;

                ptarray_insert_point(pas[j], &p, pas[j]->npoints);
            }
        }

        lwpoly = lwpoly_construct(srid, NULL, multi->poly[i].rings, pas);
        lwmpoly_add_lwpoly(lwmpoly, lwpoly);
    }

    result = gserialized_from_lwgeom(lwmpoly_as_lwgeom(lwmpoly), &ret_size);
    SET_VARSIZE(result, ret_size);

    PG_RETURN_POINTER(result);
}

PG_FUNCTION_INFO_V1(concave_hull_finalfn);
Datum concave_hull_finalfn(PG_FUNCTION_ARGS) {
    Datum geometry_array = 0;
    Datum result = 0;
    int dims[1];
    int lbs[1];
    pgis_abs *p;

    if (PG_ARGISNULL(0))
        PG_RETURN_NULL();

    p = (pgis_abs *) PG_GETARG_POINTER(0);

    dims[0] = p->a->nelems;
    lbs[0] = 1;
    geometry_array = makeMdArrayResult(p->a, 1, dims, lbs, CurrentMemoryContext, false);

    result = DirectFunctionCall1(concave_hull, geometry_array);

    if (! result)
        PG_RETURN_NULL();

    PG_RETURN_DATUM(result);
}
