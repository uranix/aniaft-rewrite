#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "struct2.h"
#include "tree2.h"
#include "memory2.h"
#include "user2.h"
#include "region2.h"

StrucMesh2  mesh;
extern StrucTree2  tree2;
extern int boolFAF;
int boolSize=0,boolNav=0,boolRegul=0;

extern int boolFAFglobal, meshnSmoothglobal;

void initParameter( void)
{
    char  *dummy=0;

    boolFAF = boolFAFglobal;
    boolSize = 1;
    mesh.nSmooth = meshnSmoothglobal;
    initUserParameters(dummy);
    return;
}  /*initParameter*/


void init( void )
{
    /* printf(" Generator  of  unstructured  2-D  meshes.\n");
       printf(" (c) Copyright  1993-95.\n");
       printf(" LabNumMath ( INM ).\n");
       printf(" Version July, 18, 1995.\n");
       */
    mesh.nPoint = 0;
    mesh.nRegion = 1;
    mesh.nRLine = 0;
    mesh.iRLine = 0;
    mesh.nSmooth = 5;
    strcpy(mesh.outFileName,"");
    strcpy(mesh.inFileName,"");
    strcpy(mesh.debFileName,"debug");

    mesh.debug = 1;
    tree2.nFace = 0;

    initMemory();
    initParameter();

    tree2.root->entrycount = 0;
    tree2.root->parent = 0;
    tree2.root->nodelist[0] = 0;
    tree2.root->nodelist[1] = 0;
    tree2.root->nodelist[2] = 0;
    tree2.root->nodelist[3] = 0;
    tree2.root->firstentry = 0;
    return;
}

void addPoint(double x, double y) {
    mesh.pts.push_back(Point(x, y));
}

void addTria(int v1, int v2, int v3, int lab) {
    mesh.tri.push_back(Triangle{v1, v2, v3, lab});
}

extern int    nVRTglobal;
extern double *vrtglobal;
extern int    nTRIglobal, *triglobal, *labtriglobal;

void outMesh() {
    nVRTglobal = mesh.nPoint;
    vrtglobal = (double*) malloc(sizeof(double)*(2*nVRTglobal));

    for (int i=0; i<mesh.nPoint; i++) {
        vrtglobal[2*i]   = mesh.pts[i].x;
        vrtglobal[2*i+1] = mesh.pts[i].y;
    }

    nTRIglobal = mesh.tri.size();
    triglobal = (int*) malloc(sizeof(int)*(3*nTRIglobal));
    labtriglobal = (int*) malloc(sizeof(int)*(nTRIglobal));
    for(size_t i=0; i<mesh.tri.size(); i++) {
        triglobal[3*i]   = mesh.tri[i].v1 + 1;
        triglobal[3*i+1] = mesh.tri[i].v2 + 1;
        triglobal[3*i+2] = mesh.tri[i].v3 + 1;
        labtriglobal[i] = mesh.tri[i].label;
    }

    if( boolNav )
        userOutMesh();
}
