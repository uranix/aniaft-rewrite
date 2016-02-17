#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#include"memory.h"
#include"tree.h"


/* extern  variables */
extern StrucMesh2  mesh;
extern int boolFAF;

/* global  variables */
StrucTree2  tree;


/* static  function */
static void addFaceArray( PStrucFace2  face );
static void remFaceArray( int  face );
static int  sectQuad( void );
static void vicinityFacesRec( PStrucNode2d  node );

int direction(double x, double y, double xc, double yc) {
    int dx=0, dy=0;
    if (x > xc)
        dx = 1;
    if (y > yc)
        dy = 2;
    return dx + dy;
}

void center(double *x, double *y, double side, int dir) {
    int dx=0, dy=0;

    if ((dir<0)||(dir>3))
        dir = 0;

    dx = dir%2;
    dy = dir/2;

    if (dx)
        x[0] += side;
    else
        x[0] -= side;

    if (dy)
        y[0] += side;
    else
        y[0] -= side;

    return;
}

double distance(double x,double y,double xc,double yc) {
    return sqrt(distanceS(x, y, xc, yc));
}/*distance*/

double distanceS(double x,double y,double xc,double yc) {
    return (x-xc)*(x-xc)+(y-yc)*(y-yc);
}/*distanceS*/

static void addFaceArray( PStrucFace2  face )
{
    int son,fath,i;
    PStrucFace2  cface,*tface=tree.face;
    if( tree.nFace >= tree.maxFace ) {
        fprintf(stderr, "aniAFT: no memory for new face in front\n");
        return;
    }

    tface[tree.nFace] = face;
    face->f = tree.nFace;
    son=++tree.nFace;
aaa:
    fath=son/2;

    if ((fath>0) && (tface[son-1]->s < tface[fath-1]->s)) {
        cface=tface[son-1];
        tface[son-1]=tface[fath-1];
        tface[fath-1]=cface;
        i=cface->f;  cface->f=tface[son-1]->f;  tface[son-1]->f=i;
        son=fath;
        if(son>1)
            goto aaa;
    }

    return;
} /*addFaceArray*/

static void remFaceArray( int  fath )
{
    int son1,son2,i;
    PStrucFace2 face,*tface;
    tface = tree.face;
    if( (fath >= tree.nFace) || (fath<0) )  {
        fprintf(stderr, "aniAFT: delete face: tree is empty\n");
        return;
    }

    tree.nFace--;
    tface[fath]=tface[tree.nFace];
    tface[fath]->f=fath;
aaa:i = 0;
    if (2*fath+1 < (tree.nFace-1)) {
        son1=2*fath+1;
        son2=son1++;
        if ((tface[fath]->s<=tface[son1]->s)&&(tface[fath]->s<=tface[son2]->s))
            i=0;
        if ((tface[fath]->s>=tface[son1]->s)&&(tface[son1]->s>=tface[son2]->s))
            i=son2;
        if ((tface[fath]->s>=tface[son2]->s)&&(tface[son2]->s>tface[son1]->s))
            i=son1;
        if ((tface[son1]->s>=tface[fath]->s)&&(tface[fath]->s>tface[son2]->s))
            i=son2;
        if ((tface[son2]->s>=tface[fath]->s)&&(tface[fath]->s>tface[son1]->s))
            i=son1;
    }
    if( 2*fath+1 == (tree.nFace-1) ){
        son1=2*fath+1;
        if((tface[fath]->s>tface[son1]->s)) i=son1;else i=0;
    }

    if(i)
    {  face=tface[i];  tface[i]=tface[fath];   tface[fath]=face;
        son1=face->f;  face->f=tface[i]->f;  tface[i]->f=son1;
        fath=i;   goto aaa;
    }

    return;
} /*remFaceArray*/


PStrucFace2 addFace(int v1, int v2, int twin )
{
    int         i,d;
    double      x1=-1.0,y1=-1.0;
    double      xc=0.5,yc=0.5,side=0.5;
    PStrucFace2 face;
    PStrucNode2d old_node,node;
    entrychain  *faceentry, *new_faceentry;

    node=tree.root;
    face = (StrucFace2*)nearAlloc( S_StrucFace2 );
    face->v1 = v1;
    face->v2 = v2;

    const Point &p1 = mesh.pts[v1];
    const Point &p2 = mesh.pts[v2];

    double x, y;

    if( twin == 0 ){
        x = (p1.x + p2.x)/2.;
        y = (p1.y + p2.y)/2.;
    } else {
        x = 0.49*p1.x + 0.51*p2.x;
        y = 0.49*p1.y + 0.51*p2.y;
    }
    face->x=x;
    face->y=y;
    face->s = distance(p1.x, p1.y, p2.x, p2.y);

    if ((x<0.)||(x>1.)||(y<0.)||(y>1)) {
        printf( "x= %lf,  y= %lf ", x, y );
        fprintf(stderr, "aniAFT: out of bouding box (wrong orientation in input?)\n");
    }
    while (side>face->s) {
        d = direction(x,y,xc,yc);
        side *= 0.5;
        center(&xc,&yc,side,d);
        old_node = node;
        node = node->nodelist[d];
        if (!node) {
            node = (StrucNode2d *)nearAlloc( S_StrucNode2d );
            for (i=0;i<4;i++)
                node->nodelist[i] = 0;
            node->parent = old_node;
            node->entrycount = 0;
            node->firstentry = 0;
            old_node->nodelist[d] = node;
        }
    }
    faceentry = node->firstentry;
    for (i=0;i<node->entrycount;i++) {
        x1 = faceentry->face->x;
        y1 = faceentry->face->y;
        if ( (x1==x)&&(y1==y)&&(v1==faceentry->face->v1)&&(v2==faceentry->face->v2) ) {
            printf("x,  y:  %lf, %lf\nx1, y1: %lf, %lf\nv1:   %d, v2:   %d\n v1_f: %d, v2_f: %d\n", x, y, x1, y1, v1, v2, faceentry->face->v1, faceentry->face->v2);
            fprintf(stderr, "aniAFT: dups in front (corrupted front?)\n");
        }
        faceentry = faceentry->next;
    }
    faceentry = node->firstentry;
    new_faceentry = (entrychain *)nearAlloc( S_entrychain );
    node->firstentry = new_faceentry;
    new_faceentry->face = face;
    new_faceentry->next = faceentry;
    node->entrycount++;
    addFaceArray( face );

    return face;
} /*addFace*/

void remFace( PStrucFace2  face )
{
    int         i,d;
    double      x,y;
    double      xc=0.5,yc=0.5,side=0.5;
    PStrucNode2d node, old_node;
    entrychain  *faceentry, *old_faceentry;

    remFaceArray(face->f);
    x = face->x;
    y = face->y;
    node = tree.root;
    while (side>face->s) {
        d = direction(x,y,xc,yc);
        side *= 0.5;
        center(&xc,&yc,side,d);
        node = node->nodelist[d];
        if (!node) {
            fprintf(stderr, "aniAFT: delete face: no such element\n");
            return;
        }
    }
    faceentry = node->firstentry;
    if (faceentry->face == face) {
        node->firstentry = faceentry->next;
        free(face);
        free(faceentry);
        node->entrycount--;
    } else {
        do {
            old_faceentry = faceentry;
            faceentry = faceentry->next;
        } while ( (faceentry)&&(faceentry->face != face) );
        if (faceentry) {
            old_faceentry->next = faceentry->next;
            free(face);
            free(faceentry);
            node->entrycount--;
        } else {
            fprintf(stderr, "aniAFT: delete face: no such element\n");
            return;
        }
    }
    while ( (node->entrycount==0)&&(!node->nodelist[0]&&!node->nodelist[1]&&
                !node->nodelist[2]&&!node->nodelist[3])&&(node!=tree.root) ) {
        old_node = node;
        node = node->parent;
        for (i=0;i<4;i++)
            if (node->nodelist[i]==old_node) break;
        node->nodelist[i] = 0;
        free(old_node);
    }
    return;
} /* remFace */

double nearest2( int *vert, double x, double y, double size )
{
    int          i,j,vn=0,v[2];
    double       p=0,dist=2.0;
    PStrucFace2  face;

    vicinityFaces(x,y,size);
    for(i=0;i<tree.nVicinityFace;i++){
        face = tree.vicinityFace[i];
        v[0] = face->v1;
        v[1] = face->v2;
        for(j=0;j<2;j++){
            p = distance(mesh.pts[v[j]].x,mesh.pts[v[j]].y,x,y);
            if(p < dist){
                dist = p;
                vn = v[j];
            }
        }/* for  2  vert  of  face */
    }/* for  i < tree.nVicinityFace */
    vert[0] = vn;

    return( dist );
} /*nearest2*/


static int sectQuad( void )
    /*  need  Quad     : tree.xc;  tree.yc;
        tree.side;
        Vicinity : tree.xVicinity;  tree.yVicinity;
        tree.sVicinity;
        */
{
    if( tree.xc+3*tree.side < tree.xVicinity-tree.sVicinity )  return(0);
    if( tree.xc-3*tree.side > tree.xVicinity+tree.sVicinity )  return(0);
    if( tree.yc+3*tree.side < tree.yVicinity-tree.sVicinity )  return(0);
    if( tree.yc-3*tree.side > tree.yVicinity+tree.sVicinity )  return(0);

    return(1);

} /*sectQuad*/


static void vicinityFacesRec( PStrucNode2d  node )
    /*  NEED    tree.xc=0.5;tree.yc=0.5;tree.side=0.5;
        tree.nVicinityFace=0;
        tree.sVicinity , tree.xVicinity , tree.yVicinity ,
        */
{
    int         i;
    double      x,y,s;
    entrychain  *faceentry;

    if (!node) return;

    for (faceentry = node->firstentry;faceentry;faceentry=faceentry->next) {
        x = faceentry->face->x;
        y = faceentry->face->y;
        s = distance(x,y,tree.xVicinity,tree.yVicinity)-0.5*faceentry->face->s;
        if (s<=tree.sVicinity) {
            tree.vicinityFace[tree.nVicinityFace] = faceentry->face;
            tree.nVicinityFace++;
            if (tree.nVicinityFace > tree.maxVicinityFace) {
                fprintf(stderr, "aniAFT: partial result in tree search\n");
                return;
            }
        }
    }
    tree.side *= 0.5;
    x = tree.xc;
    y = tree.yc;
    s = tree.side;
    if (node->nodelist[0]) {
        center(&tree.xc, &tree.yc, tree.side, 0);
        if (sectQuad())
            vicinityFacesRec( node->nodelist[0] );
    }
    for (i=1;i<4;i++) {
        if (!node->nodelist[i]) continue;
        tree.side = s;
        tree.xc = x;
        tree.yc = y;
        center(&tree.xc, &tree.yc, tree.side, i);
        if (sectQuad())
            vicinityFacesRec( node->nodelist[i] );
    }
} /* vicinityFacesRec */

void vicinityFaces( double x, double y, double size )
{
    tree.xc=0.5;
    tree.yc=0.5;
    tree.side=0.5;

    tree.nVicinityFace = 0;
    tree.sVicinity = size;
    tree.xVicinity = x;
    tree.yVicinity = y;

    vicinityFacesRec( tree.root );
} /* vicinityFaces */

