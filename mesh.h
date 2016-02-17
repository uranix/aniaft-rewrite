#ifndef H_MESH_MESH2D
#define H_MESH_MESH2D

#include "region.h"
#include "struct.h"

class Mesh {
    Region &reg;

    int     nSmooth;
    int     nRegion,nRLine,iRLine; /**/
    int     **boundVert;
    int     *nVert,*region,*bCut,*nRPoint,*nRTria;

    std::vector<Point> pts;
    std::vector<Triangle> tri;

    std::vector<int> vb;
    std::vector<int> ve;
    std::vector<int> tria1;
    std::vector<int> tria2;

    double angle( int v1, int v2, int v ) const;
    void test_quality();
    void smoothing();
    void delPoint(int n);
    void changePoint(int v1, int v2);
    void delTria(int n);
    void changeTria(int iTria, int v1, int v2, int v3);
    void sortNeigbor(int iVert);
    void calcNeigTria();
    void calcNeigbor();
    void calcEdge();
    void pack();
    void deletePoint();
    void mergePoint();
    void splitPoint();
    void swapEdge();
    void regularity(bool Regul = false);
public:
    const bool FAF;
    Mesh(Region &reg, const bool FAF) : reg(reg), FAF(FAF) { init(); }
    /* exported  function  function */
    void init( void );
    void addPoint( double x, double y );
    void addTria( int v1, int v2, int v3, int lab );
    void outMesh( void );
};
#endif
