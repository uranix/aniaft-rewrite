/* main_boundary_wing.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
   on Microsoft Windows system, link with libf2c.lib;
   on Linux or Unix systems, link with .../path/to/libf2c.a -lm
   or, if you install libf2c.a in a standard place, with -lf2c -lm
   -- in that order, at the end of the command line, as in
   cc *.o -lf2c -lm
   Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;

/* ====================================================================== */
/* Example of using mesh generator with an analytic geometry description */
/* ====================================================================== */
/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static integer nbv = 7;
    static integer nbl = 8;
    static doublereal bv[14]	/* was [2][7] */ = {
        0.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
        1.0, 0.0,
        0.4, 0.5,
        0.6, 0.5,
        1.0, 0.5 };
    static integer bl[56]	/* was [7][8] */ = {
    //  V v1
    //    V v2
    //      V type
    //        V dummy
    //           V label
    //              V idx left
    //                V idx right
        1,2,0,-1,-1,1,0,
        4,1,0,-1,-1,1,0,
        2,3,0,-1, 1,1,0,
        7,4,0,-1, 1,1,0,
        3,7,0,-1, 1,1,0,
        6,7,2, 0,11,1,1,
        6,5,1,-1, 2,1,2,
        5,6,1,-1, 2,1,2 };
    static doublereal bltail[16]	/* was [2][8] */ = { 0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,1.,0.,.5,.5,1. };

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
            e_wsle(void);

    /* Local variables */
    static doublereal h__;
    static integer nb, nc, nt, nv;
    extern /* Subroutine */ int graph_demo__(integer *, doublereal *, integer 
            *, integer *, char *, char *, ftnlen, ftnlen);
    static integer bnd[20000]	/* was [2][10000] */;
    static doublereal crv[20000]	/* was [2][10000] */;
    static integer tri[900000]	/* was [3][300000] */;
    static doublereal vrt[300000]	/* was [2][150000] */;
    static integer ifnc[10000], ierr;
    extern /* Subroutine */ int userboundary_();
    extern /* Subroutine */ double usersize_();
    extern integer aft2dboundary_(integer *, doublereal *, integer *, integer 
            *, doublereal *, doublereal *, integer *, doublereal *, integer *,
            integer *, integer *, integer *, integer *, integer *, integer *,
            doublereal *, integer *);
    static integer labelb[10000], labelt[300000];
    extern /* Subroutine */ int registeruserfn_(U_fp);
    extern /* Subroutine */ int registersizefn_(U_fp);

    /* Fortran I/O blocks */
    static cilist io___19 = { 0, 6, 0, 0, 0 };


    /*     nvmax   - maximum number of mesh nodes */
    /*     ntmax   - maximum number of mesh triangles */
    /*     nbmax   - maximum number of boundary edges */
    /* mesh generator data specifying domain analytically */
    /* rectangle */
    /*     double precision bv(2,4),bltail(2,4) */
    /*     integer          Nbv,Nbl,bl(7,4) */
    /*     data             Nbv/4/,Nbl/4/ */
    /*     data             bv/0.0,0.25, 0.0,0.5, 1.0,0.5, 1.0,0.25/  ! boundary nodes */
    /*     data             bl/1,2,0,0,1,1,0, 3,4,0,0,2,1,0,          ! outer boundary edges */
    /*    &                    2,3,0,0,3,1,0, 4,1,0,0,4,1,0/          ! outer boundary edges */
    /*     data             bltail/0,0, 0,0, 0,0, 0,0/                ! curved data for each outer boundary edge 
    */
    /* complement of a wing NACA0012 to the unit square */
    /* AFT2D library function */
    /* the name of the user written function (see file crv_model.c) */
    /* ====================================================================== */
    registeruserfn_((U_fp)userboundary_);
    registersizefn_((U_fp)usersize_);
    /* register the name in the libr */
    h__ = .01f;
    /* Generate quasiuniform mesh with meshstep h */
    /* mesh step of the quasi-uniform mesh */
    ierr = aft2dboundary_(&nbv, bv, &nbl, bl, bltail, &h__, &nv, vrt, &nt, 
            tri, labelt, &nb, bnd, labelb, &nc, crv, ifnc);
    /* geometric da */
    /* mesh data on */
    if (ierr != 0) {
        s_stop(" error in function aft2dboundary", (ftnlen)32);
    }
    s_wsle(&io___19);
    do_lio(&c__9, &c__1, "mesh: number of triangles/vertices ", (ftnlen)35);
    do_lio(&c__3, &c__1, (char *)&nt, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nv, (ftnlen)sizeof(integer));
    e_wsle();
    /*  The name must terminate with .ps */
    /*  Demo graphics has been activated */
    /*     call graph(nv,vrt, nt,tri, 'mesh_final.ps') */
    graph_demo__(&nv, vrt, &nt, tri, "mesh_final.ps", "Mesh built from bound"
            "ary given analytically", (ftnlen)13, (ftnlen)43);
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Main program alias */ int main_ () { MAIN__ (); return 0; }
