/*****************************************************************************/
/*                                                                           */
/*  Routines for Arbitrary Precision Floating-point Arithmetic               */
/*  and Fast Robust Geometric Predicates                                     */
/*  (predicates.c)                                                           */
/*                                                                           */
/*  May 18, 1996                                                             */
/*                                                                           */
/*  Placed in the public domain by                                           */
/*  Jonathan Richard Shewchuk                                                */
/*  School of Computer Science                                               */
/*  Carnegie Mellon University                                               */
/*  5000 Forbes Avenue                                                       */
/*  Pittsburgh, Pennsylvania  15213-3891                                     */
/*  jrs@cs.cmu.edu                                                           */
/*                                                                           */
/*  This file contains C implementation of algorithms for exact addition     */
/*    and multiplication of floating-point numbers, and predicates for       */
/*    robustly performing the orientation and incircle tests used in         */
/*    computational geometry.  The algorithms and underlying theory are      */
/*    described in Jonathan Richard Shewchuk.  "Adaptive Precision Floating- */
/*    Point Arithmetic and Fast Robust Geometric Predicates."  Technical     */
/*    Report CMU-CS-96-140, School of Computer Science, Carnegie Mellon      */
/*    University, Pittsburgh, Pennsylvania, May 1996.  (Submitted to         */
/*    Discrete & Computational Geometry.)                                    */
/*                                                                           */
/*  This file, the paper listed above, and other information are available   */
/*    from the Web page http://www.cs.cmu.edu/~quake/robust.html .           */
/*                                                                           */
/*****************************************************************************/

#ifndef _predicates_h
#define _predicates_h

/* On some machines, the exact arithmetic routines might be defeated by the  */
/*   use of internal extended precision floating-point registers.  Sometimes */
/*   this problem can be fixed by defining certain values to be volatile,    */
/*   thus forcing them to be stored to memory and rounded off.  This isn't   */
/*   a great solution, though, as it slows the arithmetic down.              */
/*                                                                           */
/* To try this out, write "#define INEXACT volatile" below.  Normally,       */
/*   however, INEXACT should be defined to be nothing.  ("#define INEXACT".) */

#define INEXACT                          /* Nothing */
/* #define INEXACT volatile */

#define REAL double                      /* float or double */
#define REALPRINT doubleprint
#define REALRAND doublerand
#define NARROWRAND narrowdoublerand
#define UNIFORMRAND uniformdoublerand

/* Which of the following two methods of finding the absolute values is      */
/*   fastest is compiler-dependent.  A few compilers can inline and optimize */
/*   the fabs() call; but most will incur the overhead of a function call,   */
/*   which is disastrously slow.  A faster way on IEEE machines might be to  */
/*   mask the appropriate bit, but that's difficult to do in C.              */

#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a))
/* #define Absolute(a)  fabs(a) */

/* Many of the operations are broken up into two pieces, a main part that    */
/*   performs an approximate operation, and a "tail" that computes the       */
/*   roundoff error of that operation.                                       */
/*                                                                           */
/* The operations Fast_Two_Sum(), Fast_Two_Diff(), Two_Sum(), Two_Diff(),    */
/*   Split(), and Two_Product() are all implemented as described in the      */
/*   reference.  Each of these macros requires certain variables to be       */
/*   defined in the calling routine.  The variables `bvirt', `c', `abig',    */
/*   `_i', `_j', `_k', `_l', `_m', and `_n' are declared `INEXACT' because   */
/*   they store the result of an operation that may incur roundoff error.    */
/*   The input parameter `x' (or the highest numbered `x_' parameter) must   */
/*   also be declared `INEXACT'.                                             */

double doublerand();
double narrowdoublerand();
double uniformdoublerand();
float floatrand();
float narrowfloatrand();
float uniformfloatrand();
void exactinit();
int grow_expansion(elen, e, b, h);
int grow_expansion_zeroelim(elen, e, b, h);
int expansion_sum(elen, e, flen, f, h);
int expansion_sum_zeroelim1(elen, e, flen, f, h);
int expansion_sum_zeroelim2(elen, e, flen, f, h);
int fast_expansion_sum(elen, e, flen, f, h);
int fast_expansion_sum_zeroelim(elen, e, flen, f, h);
int linear_expansion_sum(elen, e, flen, f, h);
int linear_expansion_sum_zeroelim(elen, e, flen, f, h);
int scale_expansion(elen, e, b, h);
int scale_expansion_zeroelim(elen, e, b, h);
int compress(elen, e, h);
REAL estimate(elen, e);
REAL orient2dfast(pa, pb, pc);
REAL orient2dexact(pa, pb, pc);
REAL orient2dslow(pa, pb, pc);
REAL orient2dadapt(pa, pb, pc, detsum);
REAL orient2d(pa, pb, pc);
REAL orient3dfast(pa, pb, pc, pd);
REAL orient3dexact(pa, pb, pc, pd);
REAL orient3dslow(pa, pb, pc, pd);
REAL orient3dadapt(pa, pb, pc, pd, permanent);
REAL orient3d(pa, pb, pc, pd);
REAL incirclefast(pa, pb, pc, pd);
REAL incircleexact(pa, pb, pc, pd);
REAL incircleslow(pa, pb, pc, pd);
REAL incircleadapt(pa, pb, pc, pd, permanent);
REAL incircle(pa, pb, pc, pd);
REAL inspherefast(pa, pb, pc, pd, pe);
REAL insphereexact(pa, pb, pc, pd, pe);
REAL insphereslow(pa, pb, pc, pd, pe);
REAL insphereadapt(pa, pb, pc, pd, pe, permanent);
REAL insphere(pa, pb, pc, pd, pe);

#endif
