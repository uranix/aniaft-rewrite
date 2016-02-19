#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "struct.h"
#include "tree.h"
#include "mesh.h"
#include "memory.h"
#include "user.h"
#include "region.h"

Mesh  mesh;
Tree  tree;

extern int boolFAFglobal, mesh2nSmoothglobal;

Mesh::Mesh() : reg(*(Region *)nullptr), x(pts), y(pts), FAF(boolFAFglobal) {
    nRegion = 1;
    nRLine = 0;
    iRLine = 0;

    nSmooth = mesh2nSmoothglobal;
}

void Mesh::addPoint(double x, double y, bool skip_neib) {
    pts.push_back(Point(x, y, skip_neib));
}

void Mesh::addTria(int v1, int v2, int v3, int lab) {
    tri.push_back(Triangle(v1, v2, v3, lab));
}

extern int    nVRTglobal;
extern double *vrtglobal;
extern int    nTRIglobal, *triglobal, *labtriglobal;

void Mesh::outMesh() const {
    nVRTglobal = pts.size();
    vrtglobal = (double*) malloc(sizeof(double)*(2*nVRTglobal));

    for (size_t i=0; i < pts.size(); i++) {
        vrtglobal[2*i]   = pts[i].x;
        vrtglobal[2*i+1] = pts[i].y;
    }

    nTRIglobal = tri.size();
    triglobal = (int*) malloc(sizeof(int)*(3*nTRIglobal));
    labtriglobal = (int*) malloc(sizeof(int)*(nTRIglobal));
    for(size_t i=0; i<tri.size(); i++) {
        triglobal[3*i]   = tri[i].v1 + 1;
        triglobal[3*i+1] = tri[i].v2 + 1;
        triglobal[3*i+2] = tri[i].v3 + 1;
        labtriglobal[i] = tri[i].label;
    }
}
