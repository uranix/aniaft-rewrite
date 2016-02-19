#include <stdlib.h>
#include <math.h>

#include "region.h"
#include "user.h"
#include "tree.h"

#include "mesh.h"

extern  Mesh        mesh;
extern  Tree        tree;

typedef void (*userfn_t) (int *, double *, double *, double *);
extern userfn_t userfn;

typedef double (*sizefn_t) (double*);
extern sizefn_t sizefn;

extern double ReferenceCrd[2];
extern double ScalingFactor;

extern int    nBNDglobal, *bndglobal;
extern int    nCRVglobal, *iFNCglobal;
extern double *crvglobal;

extern int    nVRTglobal;
extern double *vrtglobal;

int     nLine;

void makeLineAdvancedFront(int *boundVert, double *par_t, int nVert, int bNum,
        int bCond, int bCurve, int region,int bCut );
void makeAddAF(int *boundVert, int nVert, int region, int bCut, int num);

static  double  xBegin,xEnd,yBegin,yEnd; /* for  boundary0 */
static  double  tBegin,tEnd; /* for  boundary 0-4 */

#define MAX_POINT 50000

static  int     bnd[5][MAX_POINT],crv1[MAX_POINT];
static  double  crv[2][MAX_POINT];
static  int     nBnd=0,nCrv=0;


void boundary0(double t, double *x, double *y)
{
    x[0] = xBegin + (xEnd-xBegin)*t;
    y[0] = yBegin + (yEnd-yBegin)*t;

    return;
}/*boundary0*/


void boundary(int n, double t, double *x, double *y)
{
    if(n == 0) {
        boundary0(t,x,y);
    } 
    else {
        (* userfn) (&n, &t, x, y);
        x[0] = (x[0] - ReferenceCrd[0]) * ScalingFactor;
        y[0] = (y[0] - ReferenceCrd[1]) * ScalingFactor;
    }
    return;
}  // boundary


double nextT( int  b, double t )
{
    int      i=0;
    double   t0,t1,x0,y0,x1,y1,
             x,y,s1=0.,s,er;

    t0 = t;
    t1 = tEnd;
    boundary(b,t0,&x0,&y0);
    boundary(b,t1,&x1,&y1);

    s = sizeFace( 0.5*(x0+x1) , 0.5*(y0+y1) );
    er= 0.001*s;
    if(Tree::distance(x0,y0,x1,y1)<=s)
        return( -1. );

    while( fabs(s1-s) >= er ){
        t = 0.5*(t0+t1);
        boundary(b,t,&x,&y);
        s1 = Tree::distance(x0,y0,x,y);
        s = sizeFace( 0.5*(x0+x) , 0.5*(y0+y) );
        er = 0.01*s;
        if(s1>s) t1=t; else t0=t;
        i++;
        if( i > 100 ){
            if( fabs(s1-s) >= 10.*er )
                fprintf(stderr, "aniAFT: line: i>100 in  nextT\n");
            break;
        }
    }

    return  t;
}/* nextT */


void smoothingBoundary( int *vert, double *parameter, int iBound, int n)
    /*******************************
      n - number  of  vert  in  line  of boundary  including  two  v-vert
     *********************************/
{  int    i,j=0;
    double t,tp,tn,sp,sn,x,y;
    int    v,vp,vn;

    while( j < 5 ){
        for(i=n-2;i>0;i--){
            v  = vert[i];
            vp = vert[i-1];  tp = parameter[i-1];
            vn = vert[i+1];  tn = parameter[i+1];
            sp = sizeFace( 0.5*(mesh.pts[v].x+mesh.pts[vp].x) , 0.5*(mesh.pts[v].y+mesh.pts[vp].y) );
            sn = sizeFace( 0.5*(mesh.pts[v].x+mesh.pts[vn].x) , 0.5*(mesh.pts[v].y+mesh.pts[vn].y) );
            t = tp+(tn-tp)*sp/(sp+sn);
            parameter[i] = t;
            boundary(iBound,t,&x,&y);
            mesh.pts[v].move(x, y);
        }
        j++;
    }

    return;
}/* smoothingBoundary */



extern int    nBrglobal, nVrglobal, *brglobal;
extern double *vrbrglobal;


void initFAFRegion( void )
{
    int       i,k,v1,v2,nVert;
    double    x,y,s=0.;
    int 		 *ind;

    if( mesh.nRegion == 0 )
        mesh.nRegion = 1;
    /* init memory for mesh.nRLine */

    nVert = nVrglobal;
    ind = (int *)malloc(sizeof(int)*(nVert+1));
    k=0;
    for (i=0;i<nVert;i++){
        x = vrbrglobal[i*2];
        y = vrbrglobal[i*2+1];
        ind[i+1] = k+1;
        for (size_t j = 0; j < mesh.pts.size(); j++) {
            const Point &p = mesh.pts[j];
            if (x==p.x && y==p.y) {
                ind[i+1] = j+1;
                break;
            }
        }
        if (ind[i+1] == k+1) {
            mesh.addPoint(x, y, true);
            k++;
        }
    }
    nVert = nBrglobal;
    nBNDglobal = nBrglobal;
    bndglobal = (int*) malloc(sizeof(int)*(5*nBNDglobal));

    for (i=0;i<nVert;i++){
        v1 = ind[brglobal[i*2]];
        v2 = ind[brglobal[i*2+1]];
        bndglobal[5*i+0] = v1;
        bndglobal[5*i+1] = v2;
        bndglobal[5*i+2] = 0;
        bndglobal[5*i+3] = 1;
        bndglobal[5*i+4] = 1;
        v1--;  v2--;
        if (v1!=v2) {
            tree.addFace(mesh, v1,v2,0);
            const Point &p1 = mesh.pts[v1];
            const Point &p2 = mesh.pts[v2];
            s+=(p1.y*p2.x-p2.y*p1.x);
        }
    }/* for i */
    if (s<0) {
        printf("\nWarning: wrong orientation of boundary edges!\n\n");
    }
    free(ind);

    nCRVglobal = 0;
    crvglobal = (double*) malloc(sizeof(double)*(2*nCRVglobal));
    iFNCglobal= (int*) malloc(sizeof(int)*(nCRVglobal));
}/*initFAFRegion*/

extern double *bvglobal;
extern double *bltailglobal;
extern int *blglobal;
extern int nVertglobal, nLineglobal;

void initLineRegion( void )
{
    int       i,j;
    int       vBegin,vEnd,bNum,bCond,bCurve,region,bCut;
    int       nVert,nVVert;
    unsigned long   p;
    double    x,y,t,*par_t;
    int       *boundVert,*vVert;

    p = (unsigned long)MAX_POINT*sizeof(int);
    boundVert = (int *)malloc(p);
    if( boundVert == NULL )  perror("aniAFT");

    p = (unsigned long)MAX_POINT*sizeof(double);
    par_t = (double *)malloc(p);
    if( par_t == NULL )  perror("aniAFT");

    nVVert = nVertglobal;
    p = (unsigned long)nVVert*sizeof(int);
    vVert = (int *)malloc(p);
    if( vVert == NULL )  perror("aniAFT");

    for (i=0;i<nVVert;i++){
        x = bvglobal[i*2];
        y = bvglobal[i*2+1];
        mesh.addPoint(x, y, true);
        vVert[i] = mesh.pts.size() - 1;
    }
    nLine = nLineglobal;

    for (i=0;i<nLine;i++){
        vBegin = blglobal[i*7];
        vEnd   = blglobal[i*7+1];
        bNum   = blglobal[i*7+2];
        bCond  = blglobal[i*7+3];
        bCurve = blglobal[i*7+4];
        region = blglobal[i*7+5];
        bCut   = blglobal[i*7+6];

        vBegin--;  vEnd--;
        if( bCut > mesh.nRegion )
            mesh.nRegion = bCut;
        if( region > mesh.nRegion )
            mesh.nRegion = region;
        if( bCut > 1 || region > 1 )
            mesh.nRLine++;
        if( bNum > 0 ){
            x = bltailglobal[i*2];
            y = bltailglobal[i*2+1];
            tBegin = x;   tEnd = y;
            boundary(bNum,tBegin,&x,&y);
            mesh.pts[vVert[vBegin]].move(x, y);
            boundary(bNum,tEnd,&x,&y);
            mesh.pts[vVert[vEnd]].move(x, y);
        }
    }/* for i */
    if( mesh.nRegion == 0 )
        mesh.nRegion = 1;
    /* init memory for mesh.nRLine */
    if( mesh.nRLine > 0 ){
        p = (unsigned long)mesh.nRLine*sizeof(int *);
        mesh.boundVert = (int **)malloc(p);
        if( mesh.boundVert == NULL )  perror("aniAFT");
        p = (unsigned long)mesh.nRLine*sizeof(int);
        mesh.nVert = (int *)malloc(p);
        if( mesh.nVert == NULL )  perror("aniAFT");
        mesh.region = (int *)malloc(p);
        if( mesh.region == NULL )  perror("aniAFT");
        mesh.bCut = (int *)malloc(p);
        if( mesh.bCut == NULL )  perror("aniAFT");
    }
    p = (unsigned long)mesh.nRegion*sizeof(int);

    nLine = nLineglobal;


    for (i=0;i<nLine;i++){
        vBegin = blglobal[i*7];
        vEnd   = blglobal[i*7+1];
        bNum   = blglobal[i*7+2];
        bCond  = blglobal[i*7+3];
        bCurve = blglobal[i*7+4];
        region = blglobal[i*7+5];
        bCut   = blglobal[i*7+6];

        vBegin--;  vEnd--;
        xBegin = mesh.pts[vVert[vBegin]].x;
        xEnd   = mesh.pts[vVert[vEnd]]  .x;
        yBegin = mesh.pts[vVert[vBegin]].y;
        yEnd   = mesh.pts[vVert[vEnd]]  .y;
        if( bNum > 0 ){
            x = bltailglobal[i*2];
            y = bltailglobal[i*2+1];
            tBegin = x;   tEnd = y;
        }
        else{
            tBegin = 0.;   tEnd = 1.;
        }

        nVert = 0;
        par_t[nVert] = tBegin;
        boundVert[nVert++] = vVert[vBegin];
        t = tBegin;
        for(j=0;;j++){
            if( nVert >= MAX_POINT ) {
                fprintf(stderr, "aniAFT: line: max point\n");
                break;
            }
            t = nextT(bNum,t);
            if( t < 0. )
                break;
            boundary(bNum,t,&x,&y);
            mesh.addPoint(x, y, true);
            par_t[nVert] = t;
            boundVert[nVert++] = mesh.pts.size() - 1;
        }/* for(bool;;) */
        par_t[nVert] = tEnd;
        boundVert[nVert++] = vVert[vEnd];
        smoothingBoundary( boundVert, par_t, bNum, nVert );
        makeLineAdvancedFront(boundVert,par_t,nVert,bNum,bCond,bCurve,region,bCut);
    }/* for i */
    free(boundVert);
    free(vVert);
    free(par_t);
    return;
}/*initLineRegion*/


void initRegion( void )
{
    nBnd=0,nCrv=0;
    if( mesh.FAF )
        initFAFRegion();
    else
        initLineRegion();
}/*initRegion*/


void  initAddRegion( int num )
{
    int  i;

    for(i=0;i<mesh.nRLine;i++) {
        if (mesh.region[i] == num || mesh.bCut[i] == num)
            makeAddAF(
                    mesh.boundVert[i],
                    mesh.nVert[i],
                    mesh.region[i],
                    mesh.bCut[i],
                    num);
    }

    return;
}/*initAddRegion*/


void makeAddAF( int *boundVert, int nVert, int region, int bCut, int num )
{
    int       i=0;

    for (i=0;i<nVert-1;i++) {
        if( region == num )
            tree.addFace(mesh, boundVert[i],boundVert[i+1],bCut>0);
        if( bCut == num )
            tree.addFace(mesh, boundVert[i+1],boundVert[i],bCut>0);
    };

    return;
}

void  makeLineAdvancedFront( int *boundVert, double *par_t, int nVert,int bNum,
        int bCond, int bCurve, int region, int bCut )
{
    int       i=0;

    if( bCut > 1 || region > 1 ){
        i = nVert*sizeof(int);
        mesh.boundVert[mesh.iRLine] = (int *)malloc(i);
        if( mesh.boundVert[mesh.iRLine] == NULL )  perror("aniAFT");
        for(i=0;i<nVert;i++){
            mesh.boundVert[mesh.iRLine][i] = boundVert[i];
        }
        mesh.nVert[mesh.iRLine] = nVert;
        mesh.region[mesh.iRLine] = region;
        mesh.bCut[mesh.iRLine] = bCut;
        mesh.iRLine++;
    }

    for (i=0;i<nVert-1;i++){
        bnd[0][nBnd] = 1+boundVert[i];
        bnd[1][nBnd] = 1+boundVert[i+1];
        bnd[2][nBnd] = 0;
        if( bNum )
            bnd[2][nBnd] = nCrv+1;
        bnd[3][nBnd] = bCond;
        bnd[4][nBnd] = bCurve;
        nBnd++;
        if (bCut > 0) {
            bnd[0][nBnd] = 1+boundVert[i+1];
            bnd[1][nBnd] = 1+boundVert[i];
            bnd[2][nBnd] = 0;
            if( bNum )
                bnd[2][nBnd] = nCrv+1;
            bnd[3][nBnd] = bCond;
            bnd[4][nBnd] = bCurve;
            nBnd++;
        }

        if( bNum > 0 ){
            crv1[nCrv] = bNum;
            crv[0][nCrv] = par_t[i];
            crv[1][nCrv] = par_t[i+1];
            nCrv++;
            if (bCut > 0) {
                crv1[nCrv] = bNum;
                crv[0][nCrv] = par_t[i+1];
                crv[1][nCrv] = par_t[i];
                nCrv++;
            }
        }

        if( region == 1 )
            tree.addFace(mesh, boundVert[i],boundVert[i+1],bCut>0);
        if( bCut == 1 )
            tree.addFace(mesh, boundVert[i+1],boundVert[i],bCut>0);
    };/* for i */

    return;
}/* makeLineAdvancedFront */

double sizeFace( double x, double y )
{
    double xy[2];
    xy[0] = x/ScalingFactor + ReferenceCrd[0];
    xy[1] = y/ScalingFactor + ReferenceCrd[1];
    return  sizefn(xy)*ScalingFactor;

}/*sizeFace*/
