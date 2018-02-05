/*----------------------------------------------------------------
 *
 * src/hull.c
 *
 * Geometry builders for internal use
 *
 *----------------------------------------------------------------
 */

#include <postgres.h>
#include <utils/geo_decls.h>
#include <predicates.h>
#include <math.h>
#include <liblwgeom.h>
#include "gbuilder.h"

/*
 * Robust counter clockwise rotation algorithm. Counter clockwise if
 * ccw > 0, clockwise if ccw < 0 and collinear if ccw == 0.
 *
 * Uses the support of the robust predicates library for some special
 * cases where round off becomes an issue.
 */
static float8 _ccw(float8 ax, float8 ay, float8 bx, float8 by, float8 cx, float8 cy) {
    float8 detleft, detright, det;
    float8 detsum, errbound;

    detleft = (ax - cx) * (by - cy);
    detright = (ay - cy) * (bx - cx);
    det = detleft - detright;

    if (detleft > 0.0)
        if (detright <= 0.0)
            return det;
        else
            detsum = detleft + detright;
    else if (detleft < 0.0)
        if (detright >= 0.0)
            return det;
        else
            detsum = -detleft - detright;
    else
        return det;

    extern float8 ccwerrboundA;
    errbound = ccwerrboundA * detsum;
    if ((det >= errbound) || (-det >= errbound))
        return det;

    float8 pa[2] = {ax, ay};
    float8 pb[2] = {bx, by};
    float8 pc[2] = {cx, cy};

    return orient2dadapt(pa, pb, pc, detsum);
}

/*
 * Checks if point falls within the circumcircle of a triangle. The normal
 * algorithm is easy and only needs to compare the distance between the
 * center of the circumcircle and the point with the radius of the circumcircle.
 *
 * However, due to the existance of supertriangles, some extremely big
 * circumcircles can exists for some inputs. This is fixed by using the
 * the robust predicates library on anything that approaches a possible
 * containment.
 */
static bool _within(POINT2D *points, tri_t *t, int i, int k, float8 *dx, float8 *dy, float8 *r) {
    *dx = fabs(points[i].x - t->tri[k].x);
    *dy = fabs(points[i].y - t->tri[k].y);
    *r = t->tri[k].r;

    /*
     * Only triangles within the current x-axis range are checked, so
     * the first check only needs a y-axis comparison. Due to rounding
     * errors anything that might possibly fall within the circumcircle
     * is handled by the robust predicates library.
     */
    if ((*dy * *dy) > (*r + EPSILON))
        return false;

    float8 adx, bdx, cdx, ady, bdy, cdy;
    float8 bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
    float8 alift, blift, clift;
    float8 det;
    float8 permanent, errbound;

    adx = points[t->tri[k].p3].x - points[i].x;
    bdx = points[t->tri[k].p2].x - points[i].x;
    cdx = points[t->tri[k].p1].x - points[i].x;
    ady = points[t->tri[k].p3].y - points[i].y;
    bdy = points[t->tri[k].p2].y - points[i].y;
    cdy = points[t->tri[k].p1].y - points[i].y;

    bdxcdy = bdx * cdy;
    cdxbdy = cdx * bdy;
    alift = adx * adx + ady * ady;

    cdxady = cdx * ady;
    adxcdy = adx * cdy;
    blift = bdx * bdx + bdy * bdy;

    adxbdy = adx * bdy;
    bdxady = bdx * ady;
    clift = cdx * cdx + cdy * cdy;

    det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);
    permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * alift + (Absolute(cdxady) + Absolute(adxcdy)) * blift + (Absolute(adxbdy) + Absolute(bdxady)) * clift;

    extern float8 iccerrboundA;
    errbound = iccerrboundA * permanent;

    if ((det > errbound) || (-det > errbound))
        return (det >= 0);

    float8 pa[2] = {points[t->tri[k].p3].x, points[t->tri[k].p3].y};
    float8 pb[2] = {points[t->tri[k].p2].x, points[t->tri[k].p2].y};
    float8 pc[2] = {points[t->tri[k].p1].x, points[t->tri[k].p1].y};
    float8 pd[2] = {points[i].x, points[i].y};

    return (incircleadapt(pa, pb, pc, pd, permanent) >= 0);
}

static int _ebuf_compare(const void *a, const void *b) {
    int a_sum = ((Edge*)a)->p1 + ((Edge*)a)->p2;
    int b_sum = ((Edge*)b)->p1 + ((Edge*)b)->p2;

    if (a_sum == b_sum) {
        if ((((Edge*)a)->p1 == ((Edge*)b)->p1) || (((Edge*)a)->p1 == ((Edge*)b)->p2)) {
            return 0;
        } else {
            return 1;
        }
    } else {
        return (a_sum < b_sum) ? -1 : (a_sum > b_sum);
    }
}

static int _edge_compare(const void *a, const void *b) {
    if ((((hull_t*)a)->x[0] == ((hull_t*)b)->x[0]) && (((hull_t*)a)->y[0] == ((hull_t*)b)->y[0])) {
        if (((hull_t*)a)->x[1] == ((hull_t*)b)->x[1]) {
            return ((hull_t*)a)->y[1] < ((hull_t*)b)->y[1] ? -1 : ((hull_t*)a)->y[1] > ((hull_t*)b)->y[1];
        } else {
            return ((hull_t*)a)->x[1] < ((hull_t*)b)->x[1] ? -1 : ((hull_t*)a)->x[1] > ((hull_t*)b)->x[1];
        }
    } else if (((hull_t*)a)->x[0] == ((hull_t*)b)->x[0]) {
        return ((hull_t*)a)->y[0] < ((hull_t*)b)->y[0] ? -1 : ((hull_t*)a)->y[0] > ((hull_t*)b)->y[0];
    } else {
        return ((hull_t*)a)->x[0] < ((hull_t*)b)->x[0] ? -1 : ((hull_t*)a)->x[0] > ((hull_t*)b)->x[0];
    }
}

static int _int_compare(const void *a, const void *b) {
    return *(int*)a < *(int*)b ? -1 : *(int*)a > *(int*)b;
}

static int _tri_compare(const void *a, const void *b) {
    return ((Triangle*)a)->r < ((Triangle*)b)->r ? -1 : ((Triangle*)a)->r > ((Triangle*)b)->r;
}

static int _pt_compare(const void *a, const void *b) {
    if (((POINT2D*)a)->x == ((POINT2D*)b)->x)
        return ((POINT2D*)a)->y < ((POINT2D*)b)->y ? -1 : ((POINT2D*)a)->y > ((POINT2D*)b)->y;
    else
        return ((POINT2D*)a)->x < ((POINT2D*)b)->x ? -1 : ((POINT2D*)a)->x > ((POINT2D*)b)->x;
}

/*
 * Inserts a triangle into an array. It also calculates
 * the circumcircle of the triangle to optimize the
 * delaunay triangulation.
 */
static void _insert_tri(tri_t *t, POINT2D *points, int p1, int p2, int p3, bool flag, int *j, int *n) {
    if (t->len + 1 >= t->size) {
        if (! t->size) {
            t->size = 4;
            t->tri = (Triangle *) palloc(sizeof(Triangle) * 4);

        } else {
            t->size *= 2;
            t->tri = (Triangle *) repalloc_huge(t->tri, sizeof(Triangle) * t->size);
        }
    }

    /*
     * Use the janitor index when possible. This keeps
     * the total size of the array as small as possible,
     * while not having to shift the array every time.
     */
    int i = t->len;
    if (*n > 0) {
        i = j[--(*n)];
    } else {
        t->len++;
    }

    /*
     * Calculate the circumcircle of the triangle using the
     * length of the three sides. By storing the circumcircle
     * a lot of unnecessary calculations are skipped.
     */
    float8 ax, ay, cx, cy, az, cz, bm, nx, ny, tx, ty;

    ax = points[p2].x - points[p1].x;
    ay = points[p2].y - points[p1].y;
    cx = points[p3].x - points[p1].x;
    cy = points[p3].y - points[p1].y;
    az = ax * (points[p1].x + points[p2].x) + ay * (points[p1].y + points[p2].y);
    cz = cx * (points[p1].x + points[p3].x) + cy * (points[p1].y + points[p3].y);
    bm = 2.0 * (ax * (points[p3].y - points[p2].y) - ay * (points[p3].x - points[p2].x));

    nx = (cy * az - ay * cz) / bm;
    ny = (ax * cz - cx * az) / bm;
    tx = nx - points[p1].x;
    ty = ny - points[p1].y;

    /* Update the values at the new position */
    t->tri[i].x = nx;
    t->tri[i].y = ny;
    t->tri[i].r = tx * tx + ty * ty;
    t->tri[i].p1 = p1;
    t->tri[i].p2 = p2;
    t->tri[i].p3 = p3;
    t->tri[i].flag = flag;
}

/*
 * Calculates the delaunay triangulation for the set of
 * points in *g. The algorithm loops over the points from
 * west to east and a great optimization is achieved by filtering
 * out old points.
 */
void triangulate(POINT2D *points, int npoints, tri_t *t, int srid) {
    Edge *ebuf = (Edge *) palloc(sizeof(Edge) * 30);
    bool flag;
    float8 dx, dy;
    float8 r;
    int *buf = palloc(sizeof(int) * 10);
    int *jtri = palloc(sizeof(int) * 10);
    int bmax = 30;
    int i, j, k;
    int jmax = 10;
    int nedge, ntri = 1;
    int nt = 0;
    int tmax = 10;

    /* Points have to be sorted */
    qsort(points, npoints, sizeof(POINT2D), _pt_compare);

    /*
     * We construct a supertriangle by drawing a triangle
     * around the bounding box of all the points. This is
     * the first triangle in the buffer and it will be
     * filtered out before the final hull construction.
     */
    if (srid != 4326) {
        float8 xmid, ymid;
        float8 xmax = -INFINITY;
        float8 xmin = INFINITY;
        float8 ymax = -INFINITY;
        float8 ymin = INFINITY;

        for (i = 0; i < npoints; i++) {
            if (points[i].x < xmin)
                xmin = points[i].x;

            if (points[i].x > xmax)
                xmax = points[i].x;

            if (points[i].y < ymin)
                ymin = points[i].y;

            if (points[i].y > ymax)
                ymax = points[i].y;
        }

        dx = xmax - xmin;
        dy = ymax - ymin;
        xmid = (xmax + xmin) / 2.0;
        ymid = (ymax + ymin) / 2.0;

        /* The supertriangle vertices are added to the point buffer */
        if (dx > dy) {
            points[npoints+0].x = xmin - dy;
            points[npoints+0].y = ymin;

            points[npoints+1].x = xmid;
            points[npoints+1].y = ymax + dx;

            points[npoints+2].x = xmax + dy;
            points[npoints+2].y = ymin;

        } else {
            points[npoints+0].x = xmin;
            points[npoints+0].y = ymax + dx;

            points[npoints+1].x = xmax + dy;
            points[npoints+1].y = ymid;

            points[npoints+2].x = xmin;
            points[npoints+2].y = ymin - dx;
        }

    } else {
        points[npoints+0].x = -360.;
        points[npoints+0].y = -90.;

        points[npoints+1].x = 0.;
        points[npoints+1].y = 540.;

        points[npoints+2].x = 360.;
        points[npoints+2].y = -90.;
    }

    /* The first triangle is the supertriangle */
    _insert_tri(t, points, npoints+0, npoints+1, npoints+2, true, jtri, &nt);
    buf[0] = 0;

    /* Loop over all the x-axis sorted points */
    for (i = 0; i < npoints; i++) {
        /* Loop over all the potential triangles */
        nedge = 0, nt = 0;
        for (j = 0; j < ntri; j++) {
            k = buf[j];

            /*
             * We add the triangle edges to the edge buffer if the point is
             * inside the current circumcircle. When a point is within the
             * circumcircle this area needs to be updated.
             */
            if (_within(points, t, i, k, &dx, &dy, &r)) {
                if (nedge + 3 >= bmax) {
                    bmax *= 2;
                    ebuf = repalloc(ebuf, sizeof(Edge) * bmax);

                    jmax *= 2;
                    jtri = repalloc(jtri, sizeof(int) * jmax);
                }

                ebuf[nedge].p1 = t->tri[k].p1;
                ebuf[nedge++].p2 = t->tri[k].p2;
                ebuf[nedge].p1 = t->tri[k].p2;
                ebuf[nedge++].p2 = t->tri[k].p3;
                ebuf[nedge].p1 = t->tri[k].p3;
                ebuf[nedge++].p2 = t->tri[k].p1;

                /* The current triangle is removed because it will be split in multiple new triangles */
                jtri[nt++] = k;

            /*
             * The triangle is removed from the buffer if its range has
             * been exceeded. No new points are ever going to be able to
             * intersect its circumcircle so there is no point in checking.
             */
            } else if ((t->tri[k].x < points[i].x) && ((dx * dx) > (r + EPSILON))) {
                buf[j--] = buf[--ntri];
            }
        }

        /*
         * Edges are sorted so that the duplicates can be filtered
         * out of the buffer. We only need unique edges for this part.
         */
        qsort(ebuf, nedge, sizeof(Edge), _ebuf_compare);
        for (j = 0; j < nedge; j++) {
            if ((j < nedge - 1) && ((ebuf[j].p1 == ebuf[j+1].p2) && (ebuf[j].p2 == ebuf[j+1].p1))) {
                j++;
                continue;
            }

            /* The buffer is updated */
            if (nt == 0) {
                if (ntri + 1 >= tmax) {
                    tmax *= 2;
                    buf = repalloc_huge(buf, sizeof(int) * tmax);
                }
                buf[ntri++] = t->len;
            }

            /* Flag triangles based on the supertriangle vertices (these get filtered out later on) */
            flag = ((ebuf[j].p1 >= npoints) || (ebuf[j].p2 >= npoints));

            /* The new triangle is saved by either overwriting an old index or by adding an extra position */
            _insert_tri(t, points, ebuf[j].p1, ebuf[j].p2, i, flag, jtri, &nt);
        }
    }

    /* The triangles are sorted based on their circumcircle radius */
    qsort(t->tri, t->len, sizeof(Triangle), _tri_compare);
}

/*
 * Returns true if an edge intersects another edge. Special
 * cases such as collinearity are not checked because duplicate edges
 * are never given to this function.
 */
static bool _edge_intersection(float8 ax, float8 ay, float8 bx, float8 by, float8 cx, float8 cy, float8 dx, float8 dy) {
    float8 a = _ccw(ax, ay, bx, by, cx, cy);
    float8 b = _ccw(ax, ay, bx, by, dx, dy);
    float8 c = _ccw(cx, cy, dx, dy, ax, ay);
    float8 d = _ccw(cx, cy, dx, dy, bx, by);

    return ((((a > 0) - (a < 0)) != ((b > 0) - (b < 0))) && (((c > 0) - (c < 0)) != ((d > 0) - (d < 0))));
}

/*
 * Returns true if the point falls within the polygon. This function
 * uses a winding number algorithm so that inner rings are properly supported.
 */
static bool _point_in_polygon(polygon_t *poly, int ring, float8 x, float8 y) {
    int i = poly->index[ring], j, n = poly->index[ring+1];
    int wn = 0;

    for (i, j = i + 1; j < n; i = j++) {
        if (poly->pt[i].y <= y) {
            if (poly->pt[j].y > y)
                if (_ccw(poly->pt[i].x, poly->pt[i].y, poly->pt[j].x, poly->pt[j].y, x, y) > 0)
                    wn++;

        } else {
            if (poly->pt[j].y <= y)
                if (_ccw(poly->pt[i].x, poly->pt[i].y, poly->pt[j].x, poly->pt[j].y, x, y) < 0)
                    wn--;
        }
    }

    return wn;
}

/*
 * Wrapper around _point_in_polygon().
 */
static bool _point_in_inner_polygon(polygon_t *poly, float8 x, float8 y) {
    for (int i = 1; i < poly->rings; i++)
        if (_point_in_polygon(poly, i, x, y))
            return true;

    return false;
}

/*
 * Adding a new polygon consists of either adding a new
 * ring to an existing polygon or a new polygon entirely.
 */
static void _insert_new_polygon(multipolygon_t *multi, float8 x, float8 y, int *ptr) {
    for (int i = 0; i < multi->len; i++) {
        if (_point_in_polygon(&(multi->poly[i]), 0, x, y)) {

            /* Recursive inner rings are stored as a new polygon */
            if (_point_in_inner_polygon(&(multi->poly[i]), x, y))
                break;

            /* Dereference for readability */
            polygon_t *poly = &(multi->poly[i]);
            int j = poly->index[poly->rings];

            if (j + 1 >= poly->size) {
                poly->size *= 2;
                poly->pt = (Point *) repalloc(poly->pt, sizeof(Point) * poly->size);
            }

            /* Add the new inner ring */
            poly->pt[j].x = x;
            poly->pt[j].y = y;
            poly->rings += 1;
            poly->index = repalloc(poly->index, sizeof(int) * (poly->rings + 1));
            poly->index[poly->rings] = j + 1;
            *ptr = i;

            return;
        }
    }

    if (multi->len + 1 >= multi->size) {
        if (! multi->len) {
            multi->size = 4;
            multi->poly = (polygon_t *) palloc(sizeof(polygon_t) * 4);

        } else {
            multi->size *= 2;
            multi->poly = (polygon_t *) repalloc(multi->poly, sizeof(polygon_t) * multi->size);
        }
    }

    /* Create polygon buffer and add a new point */
    multi->poly[multi->len].pt = (Point *) palloc(sizeof(polygon_t) * 4);
    multi->poly[multi->len].index = palloc0(sizeof(int) * 2);
    multi->poly[multi->len].pt[0].x = x;
    multi->poly[multi->len].pt[0].y = y;
    multi->poly[multi->len].size = 4;
    multi->poly[multi->len].rings = 1;
    multi->poly[multi->len].index[1] = 1;
    *ptr = multi->len++;
}

/*
 * Finds the correct index based on the current
 * ring, checks for space and adds the point.
 */
static void _add_point_to_polygon(multipolygon_t *multi, float8 x, float8 y, int ptr, bool flag) {
    polygon_t *poly = &(multi->poly[ptr]);
    int i = poly->index[poly->rings];

    if (i + 1 >= poly->size) {
        poly->size *= 2;
        poly->pt = (Point *) repalloc(poly->pt, sizeof(Point) * poly->size);
    }

    /* Edge intersections are only possible in some special cases */
    if (flag) {
        float8 dx1, dy1, dx2, dy2, f, nx, ny, mx, my;
        int ii = i - 1, j, k, l;

        for (j = 0; j <= ptr; j++) {
            for (k = 0, l = 1; l < multi->poly[j].index[multi->poly[j].rings]; k = l++) {
                nx = poly->pt[ii].x + EPSILON * (x - poly->pt[ii].x);
                ny = poly->pt[ii].y + EPSILON * (y - poly->pt[ii].y);

                /*
                 * If there is a self intersection with the new point, it
                 * is pivoted outwards to a non self intersecting location.
                 */
                if (_edge_intersection(multi->poly[j].pt[k].x, multi->poly[j].pt[k].y, multi->poly[j].pt[l].x, multi->poly[j].pt[l].y, nx, ny, x, y)) {
                    if (_ccw(nx, ny, x, y, multi->poly[j].pt[k].x, multi->poly[j].pt[k].y) > 0) {
                        mx = multi->poly[j].pt[l].x + -EPSILON * (multi->poly[j].pt[k].x - multi->poly[j].pt[l].x);
                        my = multi->poly[j].pt[l].y + -EPSILON * (multi->poly[j].pt[k].y - multi->poly[j].pt[l].y);

                    } else {
                        mx = multi->poly[j].pt[k].x + -EPSILON * (multi->poly[j].pt[l].x - multi->poly[j].pt[k].x);
                        my = multi->poly[j].pt[k].y + -EPSILON * (multi->poly[j].pt[l].y - multi->poly[j].pt[k].y);
                    }

                    dx1 = x - nx;
                    dy1 = y - ny;
                    dx2 = mx - nx;
                    dy2 = my - ny;

                    f = 1 / (sqrt(dx2 * dx2 + dy2 * dy2) / sqrt(dx1 * dx1 + dy1 * dy1));

                    poly->pt[i].x = nx + f * (mx - nx);
                    poly->pt[i].y = ny + f * (my - ny);
                    poly->index[poly->rings]++;

                    return;
                }
            }
        }
    }

    poly->pt[i].x = x;
    poly->pt[i].y = y;
    poly->index[poly->rings]++;
}

/*
 * When a hull does not reach its starting point it is
 * most likely an artifact from the delaunay
 * triangulation. These edges are flushed from the buffer
 * by either removing a ring or a polygon.
 */
static void _flush_uncompleted_hull(multipolygon_t *multi, int ptr, bool *new_polygon) {
    polygon_t *a = &(multi->poly[ptr]);
    polygon_t *b = &(multi->poly[multi->len-1]);

    /* Remove the entire polygon when only the outer ring is left */
    if (a->rings == 1) {
        if (ptr != multi->len - 1) {
            int i = 0;
            a->rings = 0;
            a->index[0] = 0;

            /* Deep copy the points */
            for (i; i < b->index[b->rings]; i++)
                _add_point_to_polygon(multi, b->pt[i].x, b->pt[i].y, ptr, false);

            /* Reset the rings */
            a->rings = b->rings;
            a->index = repalloc(a->index, sizeof(int) * (a->rings + 1));

            /* Finally, deep copy the ring indices */
            for (i = 1; i < a->rings + 1; i++)
                a->index[i] = b->index[i];

        }
        pfree(b->pt);
        pfree(b->index);
        multi->len--;
    }

    /* Inner rings are simple */
    if (a->rings > 1)
        a->rings--;

    *new_polygon = true;
}

/*
 * A polygon or ring might still be invalid if it does
 * not meet certain constraints. If the polygon or ring
 * has less than 4 points it automatically returns true.
 * Collinearity of the first and second to last point also
 * invalidates the input polygon or ring.
 */
static bool _is_polygon_invalid(multipolygon_t *multi, int ptr) {
    polygon_t *poly = &(multi->poly[ptr]);
    int i = poly->index[poly->rings];
    int j = poly->index[poly->rings-1];

    /* Anything with less than 4 points is not even a polygon */
    if ((i - j) < 4)
        return true;

    /* If the second to last point is on the same x-axis as the first point,
     * we also check the next points in line for this intersection. Any other
     * points do not really need to be checked, because the triangulation will
     * never allow lines to cross.
     */
    if (poly->pt[i-1].x == poly->pt[j].x)
        if ((poly->pt[i-2].x == poly->pt[j].x) || (poly->pt[j+1].x == poly->pt[i].x))
            return true;


    /* Same check as above, but for the y-axis */
    if (poly->pt[i-1].y == poly->pt[j].y)
        if ((poly->pt[i-2].y == poly->pt[j].y) || (poly->pt[j+1].y == poly->pt[i].y))
            return true;

    return false;
}

/*
 * Builds a concave hull by removing all triangles based on
 * the radius of its circumcircle. The remaining
 * triangles are the concave hull, which is constructed from
 * all the unique edges that remain. The result is often a
 * multipolygon.
 */
void construct_concave_hull(POINT2D *points, int npoints, tri_t *t, multipolygon_t *multi, float8 mod) {
    bool invalid;
    bool new_polygon = true;
    float8 alpha;
    float8 fx[2], fy[2];
    float8 lp[2], lpn[2];
    float8 v1x, v1y, v2x, v2y, dot, mag, a1, a2;
    float8 xrange;
    float8 xmin, xmax;
    hull_t *buf;
    int a, b;
    int hmax = 0;
    int i, j = 0, k, l;
    int ns = 0;
    int p1, p2, p3;
    int ptr = 0;
    int *hull;

    /*
     * Force the use of MaxAllocHugeSize instead of MaxAllocSize by immediately
     * calling repalloc_huge. There is no such function for palloc, so we do not
     * allocate everything in one go.
     */
    buf = (hull_t *) palloc(sizeof(hull_t));
    buf = (hull_t *) repalloc_huge(buf, sizeof(hull_t) * t->len * 3);

    hull = palloc(sizeof(int) * npoints);

    /* Find the range of the input points */
    for (i = 0; i < npoints; i++) {
        if (points[i].x < xmin)
            xmin = points[i].x;

        if (points[i].x > xmax)
            xmax = points[i].x;
    }

    alpha = 0.001 * mod;
    xrange = xmax - xmin;

    /*
     * We remove all the triangles based on their alpha. Remaining
     * triangles are saved as edges, which are reorientated towards the east.
     */
    for (i = 0; i < t->len; i++) {
        p1 = t->tri[i].p1;
        p2 = t->tri[i].p2;
        p3 = t->tri[i].p3;

        if ((! t->tri[i].flag) && ((t->tri[i].r / xrange) <= alpha)) {
            if (points[p1].x < points[p2].x) {
                a = p1;
                b = p2;
            } else if ((points[p1].x == points[p2].x) && (points[p1].y < points[p2].y)) {
                a = p1;
                b = p2;
            } else {
                a = p2;
                b = p1;
            }

            buf[j].x[0] = points[a].x;
            buf[j].y[0] = points[a].y;
            buf[j].x[1] = points[b].x;
            buf[j].y[1] = points[b].y;
            buf[j++].flag = false;

            if (points[p2].x < points[p3].x) {
                a = p2;
                b = p3;
            } else if ((points[p2].x == points[p3].x) && (points[p2].y < points[p3].y)) {
                a = p2;
                b = p3;
            } else {
                a = p3;
                b = p2;
            }

            buf[j].x[0] = points[a].x;
            buf[j].y[0] = points[a].y;
            buf[j].x[1] = points[b].x;
            buf[j].y[1] = points[b].y;
            buf[j++].flag = false;

            if (points[p3].x < points[p1].x) {
                a = p3;
                b = p1;
            } else if ((points[p3].x == points[p1].x) && (points[p3].y < points[p1].y)) {
                a = p3;
                b = p1;
            } else {
                a = p1;
                b = p3;
            }

            buf[j].x[0] = points[a].x;
            buf[j].y[0] = points[a].y;
            buf[j].x[1] = points[b].x;
            buf[j].y[1] = points[b].y;
            buf[j++].flag = false;
        }
    }

    /* Save all non duplicate edges */
    qsort(buf, j, sizeof(hull_t), _edge_compare);
    for (i = 0, k = 1; i < j - 1; i++, k++) {
        if ((buf[i].x[0] == buf[k].x[0]) && (buf[i].y[0] == buf[k].y[0]) && (buf[i].x[1] == buf[k].x[1]) && (buf[i].y[1] == buf[k].y[1])) {
            i++;
            k++;
            continue;
        }

        hull[hmax++] = i;
    }

    if (i < j)
        hull[hmax++] = i;

    int seen[hmax];

    while (hmax > 0) {
        if (new_polygon) {
            /* Resort the hull buffer so that the most western point can be grabbed */
            qsort(hull, hmax, sizeof(int), _int_compare);

            /* The bottom hull is always first */
            i = (_ccw(buf[hull[0]].x[0], buf[hull[0]].y[0], buf[hull[0]].x[1], buf[hull[0]].y[1], buf[hull[1]].x[1], buf[hull[1]].y[1]) < 0);

            /* Checkpoints for the edge finder */
            fx[0] = lp[0] = buf[hull[i]].x[0];
            fy[0] = lp[1] = buf[hull[i]].y[0];
            fx[1] = buf[hull[i]].x[0] + EPSILON * (buf[hull[i]].x[1] - buf[hull[i]].x[0]);
            fy[1] = buf[hull[i]].y[0] + EPSILON * (buf[hull[i]].y[1] - buf[hull[i]].y[0]);
            lpn[0] = buf[hull[i]].x[1];
            lpn[1] = buf[hull[i]].y[1];

            /* New polygons can also be inner rings (see: _insert_new_polygon()) */
            _insert_new_polygon(multi, fx[1], fy[1], &ptr);

            hull[i] = hull[--hmax];
            new_polygon = false;
        }

        /*
         * Everything in the buffer is checked. This is unfortunately required
         * to guarantee valid polygons (no self intersections).
         */
        ns = 0;
        invalid = true;
        for (i = 0; i < hmax; i++) {
            if ((buf[hull[i]].x[0] == lpn[0]) && (buf[hull[i]].y[0] == lpn[1]))
                seen[ns++] = i;
            else if ((buf[hull[i]].x[1] == lpn[0]) && (buf[hull[i]].y[1] == lpn[1]))
                seen[ns++] = i * -1 - 1;
        }

        /*
         * Multiple candidates have been seen so the correct edge needs to be
         * found. The correct edge is always the rightmost one, because this will
         * never allow a self intersection.
         */
        if (ns > 1) {
            i = seen[--ns];
            k = 1;
            v1x = lp[0] - lpn[0];
            v1y = lp[1] - lpn[1];

            /* A negative seen index is given to flipped edges */
            if (i < 0) {
                i = (i + 1) * -1;
                k = 0;
            }

            /*
             * Any duplicate is a potential case for self intersection. The flag will
             * trigger a validation when this point is ever added to the polygon.
             */
            buf[hull[i]].flag = true;

            /* Get the angle of the two edges */
            v2x = buf[hull[i]].x[k] - lpn[0];
            v2y = buf[hull[i]].y[k] - lpn[1];
            dot = v1x * v2x + v1y * v2y;
            mag = sqrt(v1x * v1x + v1y * v1y) * sqrt(v2x * v2x + v2y * v2y);
            a1 = acos(dot / mag);

            /* Inner angles are sometimes possible, but they should not overpower outer angles */
            if (_ccw(lp[0], lp[1], lpn[0], lpn[1], buf[hull[i]].x[k], buf[hull[i]].y[k]) > 0)
                a1 = (M_PI - a1) + M_PI;

            /* Every possible candidate is compared to the rest based on the angle from the previous edge */
            while (ns > 0) {
                j = seen[--ns];
                l = 1;

                if (j < 0) {
                    j = (j + 1) * -1;
                    l = 0;
                }

                buf[hull[j]].flag = true;

                v2x = buf[hull[j]].x[l] - lpn[0];
                v2y = buf[hull[j]].y[l] - lpn[1];
                dot = v1x * v2x + v1y * v2y;
                mag = sqrt(v1x * v1x + v1y * v1y) * sqrt(v2x * v2x + v2y * v2y);
                a2 = acos(dot / mag);

                if (_ccw(lp[0], lp[1], lpn[0], lpn[1], buf[hull[j]].x[l], buf[hull[j]].y[l]) > 0)
                    a2 = (M_PI - a2) + M_PI;

                if (a2 < a1) {
                    i = j;
                    k = l;
                    a1 = a2;
                }
            }

            /* Saved point is shifted a tiny bit so that a duplicate point does not cause a self intersection */
            _add_point_to_polygon(multi, buf[hull[i]].x[!k] + EPSILON * (buf[hull[i]].x[k] - buf[hull[i]].x[!k]), buf[hull[i]].y[!k] + EPSILON * (buf[hull[i]].y[k] - buf[hull[i]].y[!k]), ptr, true);

            lp[0] = lpn[0];
            lp[1] = lpn[1];
            lpn[0] = buf[hull[i]].x[k];
            lpn[1] = buf[hull[i]].y[k];

            /* If the next point is also the first (last) point, this polygon is marked as done */
            if ((lpn[0] == fx[0]) && (lpn[1] == fy[0])) {
                _add_point_to_polygon(multi, fx[1], fy[1], ptr, false);
                new_polygon = true;
            }

            hull[i] = hull[--hmax];
            invalid = (((new_polygon) && (_is_polygon_invalid(multi, ptr))) || ((! new_polygon) && (hmax == 0)));

        /* If there is only one candiate it can immediately be added to the polygon */
        } else if (ns == 1) {
            i = seen[--ns];

            /*
             * The point is added based on which side of the edge was encountered. This point is
             * then also shifted a bit to help with self intersections.
             */
            if (i < 0) {
                i = (i + 1) * -1;
                _add_point_to_polygon(multi, buf[hull[i]].x[1] + EPSILON * (buf[hull[i]].x[0] - buf[hull[i]].x[1]), buf[hull[i]].y[1] + EPSILON * (buf[hull[i]].y[0] - buf[hull[i]].y[1]), ptr, true);

                lp[0] = lpn[0];
                lp[1] = lpn[1];
                lpn[0] = buf[hull[i]].x[0];
                lpn[1] = buf[hull[i]].y[0];

            } else {
                _add_point_to_polygon(multi, buf[hull[i]].x[0] + EPSILON * (buf[hull[i]].x[1] - buf[hull[i]].x[0]), buf[hull[i]].y[0] + EPSILON * (buf[hull[i]].y[1] - buf[hull[i]].y[0]), ptr, true);

                lp[0] = lpn[0];
                lp[1] = lpn[1];
                lpn[0] = buf[hull[i]].x[1];
                lpn[1] = buf[hull[i]].y[1];
            }

            /* If the next point is also the first (last) point, this polygon is marked as done */
            if ((lpn[0] == fx[0]) && (lpn[1] == fy[0])) {
                _add_point_to_polygon(multi, fx[1], fy[1], ptr, false);
                new_polygon = true;
            }

            hull[i] = hull[--hmax];
            invalid = (((new_polygon) && (_is_polygon_invalid(multi, ptr))) || ((! new_polygon) && (hmax == 0)));
        }

        if (invalid)
            _flush_uncompleted_hull(multi, ptr, &new_polygon);
    }
}
