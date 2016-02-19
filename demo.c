/* demo.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;
static integer c__5 = 5;
static integer c__9 = 9;

/* ====================================================================== */
/* Subroutine */ int graph_demo__(integer *nv, doublereal *vrt, integer *nt, 
        integer *tri, char *fname, char *demo_message__, ftnlen fname_len, 
        ftnlen demo_message_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
            ftnlen), e_wsfe(void), s_wsle(cilist *), do_lio(integer *, 
                integer *, char *, ftnlen), e_wsle(void), f_clos(cllist *);

    /* Local variables */
    static char fnameext[30];
    static integer i__, k, l;
    static doublereal kx, scale;
    extern /* Subroutine */ int headerps_demo__(integer *);
    static doublereal xymin[2], xymax[2];

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 1, 0, "(A)", 0 };
    static cilist io___8 = { 0, 1, 0, "(A)", 0 };
    static cilist io___9 = { 0, 1, 0, "(A)", 0 };
    static cilist io___10 = { 0, 1, 0, "(A)", 0 };
    static cilist io___11 = { 0, 1, 0, "(A)", 0 };
    static cilist io___13 = { 0, 1, 0, 0, 0 };
    static cilist io___15 = { 0, 1, 0, "(3A)", 0 };
    static cilist io___16 = { 0, 1, 0, 0, 0 };


    /* ====================================================================== */
    /*  Make a simple ps-figure of the triangulation. The file name */
    /*  MUST have extension '.ps'! */
    /* ====================================================================== */
    /* ====================================================================== */
    /* Parameter adjustments */
    tri -= 4;
    vrt -= 3;

    /* Function Body */
    i__ = 1;
    while(s_cmp(fname + (i__ - 1), ".ps", (ftnlen)3, (ftnlen)3) != 0) {
        ++i__;
    }
    s_copy(fnameext, fname, (ftnlen)30, i__ + 2);
    xymin[0] = vrt[3];
    xymin[1] = vrt[4];
    xymax[0] = vrt[3];
    xymax[1] = vrt[4];
    i__1 = *nv;
    for (i__ = 2; i__ <= i__1; ++i__) {
        /* Computing MIN */
        d__1 = xymin[0], d__2 = vrt[(i__ << 1) + 1];
        xymin[0] = min(d__1,d__2);
        /* Computing MIN */
        d__1 = xymin[1], d__2 = vrt[(i__ << 1) + 2];
        xymin[1] = min(d__1,d__2);
        /* Computing MAX */
        d__1 = xymax[0], d__2 = vrt[(i__ << 1) + 1];
        xymax[0] = max(d__1,d__2);
        /* Computing MAX */
        d__1 = xymax[1], d__2 = vrt[(i__ << 1) + 2];
        xymax[1] = max(d__1,d__2);
    }
    /* Computing MAX */
    d__1 = xymax[0] - xymin[0], d__2 = xymax[1] - xymin[1];
    scale = max(d__1,d__2);
    kx = 500.f / scale;
    o__1.oerr = 0;
    o__1.ounit = 1;
    o__1.ofnmlen = fname_len;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_wsfe(&io___7);
    do_fio(&c__1, "%!PS-Adobe-2.0 EPSF-2.0", (ftnlen)23);
    e_wsfe();
    s_wsfe(&io___8);
    do_fio(&c__1, "%%BoundingBox: 0 0  520 520", (ftnlen)27);
    e_wsfe();
    s_wsfe(&io___9);
    do_fio(&c__1, "%%EndComments", (ftnlen)13);
    e_wsfe();
    s_wsfe(&io___10);
    do_fio(&c__1, " 10 10 translate 0 setlinewidth", (ftnlen)31);
    e_wsfe();
    s_wsfe(&io___11);
    do_fio(&c__1, " /t{newpath moveto lineto lineto closepath stroke}def", (
                ftnlen)53);
    e_wsfe();
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
        s_wsle(&io___13);
        for (l = 1; l <= 3; ++l) {
            for (i__ = 1; i__ <= 2; ++i__) {
                d__1 = (vrt[i__ + (tri[l + k * 3] << 1)] - xymin[i__ - 1]) * 
                    kx;
                do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal)
                      );
            }
        }
        do_lio(&c__9, &c__1, " t", (ftnlen)2);
        e_wsle();
    }
    /* ... demo part */
    headerps_demo__(&c__1);
    s_wsfe(&io___15);
    do_fio(&c__1, "250 450 (", (ftnlen)9);
    do_fio(&c__1, demo_message__, demo_message_len);
    do_fio(&c__1, ") clearANDctext", (ftnlen)15);
    e_wsfe();
    s_wsle(&io___16);
    do_lio(&c__9, &c__1, " showpage", (ftnlen)9);
    e_wsle();
    cl__1.cerr = 0;
    cl__1.cunit = 1;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* graph_demo__ */

/* ====================================================================== */
/* Subroutine */ int isolines_demo__(doublereal *u, integer *nv, doublereal *
        vrt, integer *nt, integer *tri, integer *nb, integer *bnd, char *
        fname, integer *ni, char *demo_message__, ftnlen fname_len, ftnlen 
        demo_message_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
            ftnlen), e_wsfe(void), s_wsle(cilist *), do_lio(integer *, 
                integer *, char *, ftnlen), e_wsle(void), f_clos(cllist *);

    /* Local variables */
    static char fnameext[30];
    static integer i__, j, k;
    static doublereal x[2], y[2];
    static integer jp, kp, is;
    static doublereal kx, tet, umin, umax, ucur, scale;
    extern /* Subroutine */ int headerps_demo__(integer *);
    static doublereal xymin[2], xymax[2];
    static integer icount;

    /* Fortran I/O blocks */
    static cilist io___25 = { 0, 1, 0, "(A)", 0 };
    static cilist io___26 = { 0, 1, 0, "(A)", 0 };
    static cilist io___27 = { 0, 1, 0, "(A)", 0 };
    static cilist io___28 = { 0, 1, 0, "(A)", 0 };
    static cilist io___29 = { 0, 1, 0, "(A)", 0 };
    static cilist io___41 = { 0, 1, 0, "(4F7.1,A)", 0 };
    static cilist io___42 = { 0, 1, 0, "(4F7.1,A)", 0 };
    static cilist io___43 = { 0, 1, 0, "(3A)", 0 };
    static cilist io___44 = { 0, 1, 0, 0, 0 };


    /* ====================================================================== */
    /*  Routine draws solution isolines. The number of isolines is ni. */
    /*  The file fName must have extension '.ps'. */
    /* ================================================================ */
    /* ================================================================ */
    /* Parameter adjustments */
    bnd -= 3;
    tri -= 4;
    vrt -= 3;
    --u;

    /* Function Body */
    i__ = 1;
    while(s_cmp(fname + (i__ - 1), ".ps", (ftnlen)3, (ftnlen)3) != 0) {
        ++i__;
    }
    s_copy(fnameext, fname, (ftnlen)30, i__ + 2);
    xymin[0] = vrt[3];
    xymin[1] = vrt[4];
    xymax[0] = vrt[3];
    xymax[1] = vrt[4];
    umin = u[1];
    umax = u[1];
    i__1 = *nv;
    for (i__ = 2; i__ <= i__1; ++i__) {
        /* Computing MIN */
        d__1 = xymin[0], d__2 = vrt[(i__ << 1) + 1];
        xymin[0] = min(d__1,d__2);
        /* Computing MIN */
        d__1 = xymin[1], d__2 = vrt[(i__ << 1) + 2];
        xymin[1] = min(d__1,d__2);
        /* Computing MAX */
        d__1 = xymax[0], d__2 = vrt[(i__ << 1) + 1];
        xymax[0] = max(d__1,d__2);
        /* Computing MAX */
        d__1 = xymax[1], d__2 = vrt[(i__ << 1) + 2];
        xymax[1] = max(d__1,d__2);
        /* Computing MAX */
        d__1 = umax, d__2 = u[i__];
        umax = max(d__1,d__2);
        /* Computing MIN */
        d__1 = umin, d__2 = u[i__];
        umin = min(d__1,d__2);
    }
    /* Computing MAX */
    d__1 = xymax[0] - xymin[0], d__2 = xymax[1] - xymin[1];
    scale = max(d__1,d__2);
    kx = 500.f / scale;
    o__1.oerr = 0;
    o__1.ounit = 1;
    o__1.ofnmlen = 30;
    o__1.ofnm = fnameext;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_wsfe(&io___25);
    do_fio(&c__1, "%!PS-Adobe-2.0 EPSF-2.0", (ftnlen)23);
    e_wsfe();
    s_wsfe(&io___26);
    do_fio(&c__1, "%%BoundingBox: 0 0  520 520", (ftnlen)27);
    e_wsfe();
    s_wsfe(&io___27);
    do_fio(&c__1, "%%EndComments", (ftnlen)13);
    e_wsfe();
    s_wsfe(&io___28);
    do_fio(&c__1, " 10 10 translate 0 setlinewidth", (ftnlen)31);
    e_wsfe();
    s_wsfe(&io___29);
    do_fio(&c__1, " /v{moveto lineto stroke}def", (ftnlen)28);
    e_wsfe();
    i__1 = *nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
        i__2 = *ni;
        for (is = 1; is <= i__2; ++is) {
            icount = 0;
            ucur = umin + (is - 1) * (umax - umin) / (*ni - 1);
            for (j = 1; j <= 3; ++j) {
                k = j - 1;
                if (k == 0) {
                    k = 3;
                }
                jp = tri[j + i__ * 3];
                kp = tri[k + i__ * 3];
                if (u[jp] == u[kp] && u[jp] == ucur && icount < 2) {
                    x[0] = vrt[(jp << 1) + 1];
                    y[0] = vrt[(jp << 1) + 2];
                    x[1] = vrt[(kp << 1) + 1];
                    y[1] = vrt[(kp << 1) + 2];
                    icount = 2;
                }
                if (((u[jp] < ucur && ucur <= u[kp]) || (u[jp] >= ucur && ucur > 
                                u[kp])) && icount < 2) {
                    ++icount;
                    tet = (ucur - u[jp]) / (u[kp] - u[jp]);
                    x[icount - 1] = vrt[(jp << 1) + 1] + tet * (vrt[(kp << 1) 
                            + 1] - vrt[(jp << 1) + 1]);
                    y[icount - 1] = vrt[(jp << 1) + 2] + tet * (vrt[(kp << 1) 
                            + 2] - vrt[(jp << 1) + 2]);
                }
            }
            if (icount != 0) {
                s_wsfe(&io___41);
                d__1 = (x[0] - xymin[0]) * kx;
                do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
                d__2 = (y[0] - xymin[1]) * kx;
                do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
                d__3 = (x[1] - xymin[0]) * kx;
                do_fio(&c__1, (char *)&d__3, (ftnlen)sizeof(doublereal));
                d__4 = (y[1] - xymin[1]) * kx;
                do_fio(&c__1, (char *)&d__4, (ftnlen)sizeof(doublereal));
                do_fio(&c__1, " v", (ftnlen)2);
                e_wsfe();
            }
        }
    }
    /*    draw boundary */
    i__1 = *nb;
    for (i__ = 1; i__ <= i__1; ++i__) {
        jp = bnd[(i__ << 1) + 1];
        kp = bnd[(i__ << 1) + 2];
        s_wsfe(&io___42);
        d__1 = (vrt[(jp << 1) + 1] - xymin[0]) * kx;
        do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
        d__2 = (vrt[(jp << 1) + 2] - xymin[1]) * kx;
        do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
        d__3 = (vrt[(kp << 1) + 1] - xymin[0]) * kx;
        do_fio(&c__1, (char *)&d__3, (ftnlen)sizeof(doublereal));
        d__4 = (vrt[(kp << 1) + 2] - xymin[1]) * kx;
        do_fio(&c__1, (char *)&d__4, (ftnlen)sizeof(doublereal));
        do_fio(&c__1, " v", (ftnlen)2);
        e_wsfe();
    }
    /* ... demo part */
    headerps_demo__(&c__1);
    s_wsfe(&io___43);
    do_fio(&c__1, "250 450 (", (ftnlen)9);
    do_fio(&c__1, demo_message__, demo_message_len);
    do_fio(&c__1, ") clearANDctext", (ftnlen)15);
    e_wsfe();
    s_wsle(&io___44);
    do_lio(&c__9, &c__1, " showpage", (ftnlen)9);
    e_wsle();
    cl__1.cerr = 0;
    cl__1.cunit = 1;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* isolines_demo__ */

/* ====================================================================== */
/* Subroutine */ int headerps_demo__(integer *io)
{
    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Fortran I/O blocks */
    static cilist io___45 = { 0, 0, 0, "(A)", 0 };
    static cilist io___46 = { 0, 0, 0, "(A)", 0 };
    static cilist io___47 = { 0, 0, 0, "(A,/)", 0 };
    static cilist io___48 = { 0, 0, 0, "(A)", 0 };
    static cilist io___49 = { 0, 0, 0, "(A,/)", 0 };
    static cilist io___50 = { 0, 0, 0, "(A)", 0 };
    static cilist io___51 = { 0, 0, 0, "(A,/)", 0 };
    static cilist io___52 = { 0, 0, 0, "(A)", 0 };
    static cilist io___53 = { 0, 0, 0, "(A)", 0 };
    static cilist io___54 = { 0, 0, 0, "(A)", 0 };
    static cilist io___55 = { 0, 0, 0, "(A)", 0 };
    static cilist io___56 = { 0, 0, 0, "(A,/)", 0 };
    static cilist io___57 = { 0, 0, 0, "(A)", 0 };
    static cilist io___58 = { 0, 0, 0, "(A)", 0 };
    static cilist io___59 = { 0, 0, 0, "(A)", 0 };
    static cilist io___60 = { 0, 0, 0, "(A,/)", 0 };
    static cilist io___61 = { 0, 0, 0, "(A)", 0 };


    /* ====================================================================== */
    /* The routines writes a header for a PS demo file */
    /* ====================================================================== */
    io___45.ciunit = *io;
    s_wsfe(&io___45);
    do_fio(&c__1, "/l{lineto}def /m{moveto}def /s{l stroke}def", (ftnlen)43);
    e_wsfe();
    io___46.ciunit = *io;
    s_wsfe(&io___46);
    do_fio(&c__1, "/x{exch}def /rm{rmoveto}def", (ftnlen)27);
    e_wsfe();
    io___47.ciunit = *io;
    s_wsfe(&io___47);
    do_fio(&c__1, "/np{newpath}def /cp{closepath}def", (ftnlen)33);
    e_wsfe();
    io___48.ciunit = *io;
    s_wsfe(&io___48);
    do_fio(&c__1, "/fsize 27 def", (ftnlen)13);
    e_wsfe();
    io___49.ciunit = *io;
    s_wsfe(&io___49);
    do_fio(&c__1, "/Times-Roman findfont fsize scalefont setfont", (ftnlen)45)
        ;
    e_wsfe();
    io___50.ciunit = *io;
    s_wsfe(&io___50);
    do_fio(&c__1, "/ctext{3 1 roll m dup stringwidth fsize 0.33", (ftnlen)44);
    e_wsfe();
    io___51.ciunit = *io;
    s_wsfe(&io___51);
    do_fio(&c__1, "mul sub x 2 div 0 x sub x rm show}def", (ftnlen)37);
    e_wsfe();
    io___52.ciunit = *io;
    s_wsfe(&io___52);
    do_fio(&c__1, "/clearField{/yUp x def /xRight x def", (ftnlen)36);
    e_wsfe();
    io___53.ciunit = *io;
    s_wsfe(&io___53);
    do_fio(&c__1, "/yDown x def /xLeft x def", (ftnlen)25);
    e_wsfe();
    io___54.ciunit = *io;
    s_wsfe(&io___54);
    do_fio(&c__1, "gsave np xLeft yDown m xRight yDown", (ftnlen)35);
    e_wsfe();
    io___55.ciunit = *io;
    s_wsfe(&io___55);
    do_fio(&c__1, "l xRight yUp l xLeft yUp l cp 1 setgray fill", (ftnlen)44);
    e_wsfe();
    io___56.ciunit = *io;
    s_wsfe(&io___56);
    do_fio(&c__1, "grestore}def", (ftnlen)12);
    e_wsfe();
    io___57.ciunit = *io;
    s_wsfe(&io___57);
    do_fio(&c__1, "/clearANDctext{3 copy stringwidth /dY fsize 2", (ftnlen)45)
        ;
    e_wsfe();
    io___58.ciunit = *io;
    s_wsfe(&io___58);
    do_fio(&c__1, "div def pop /dX x 2 div def /yDown x def /xLeft", (ftnlen)
            47);
    e_wsfe();
    io___59.ciunit = *io;
    s_wsfe(&io___59);
    do_fio(&c__1, "x def xLeft dX sub yDown dY sub xLeft dX add", (ftnlen)44);
    e_wsfe();
    io___60.ciunit = *io;
    s_wsfe(&io___60);
    do_fio(&c__1, "yDown dY add clearField ctext}def", (ftnlen)33);
    e_wsfe();
    io___61.ciunit = *io;
    s_wsfe(&io___61);
    do_fio(&c__1, "0. 0. 1. setrgbcolor", (ftnlen)20);
    e_wsfe();
    return 0;
} /* headerps_demo__ */

