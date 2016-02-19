#ifndef H_MESH_MESH2D
#define H_MESH_MESH2D

#include "region.h"
#include "struct.h"

extern int boolFAFglobal;

struct XCoord {
    std::vector<Point> &pts;
    XCoord(std::vector<Point> &pts) : pts(pts) { }
    double &operator[](int i) { return pts[i].x; }
    const double &operator[](int i) const { return pts[i].x; }
};

struct YCoord {
    std::vector<Point> &pts;
    YCoord(std::vector<Point> &pts) : pts(pts) { }
    double &operator[](int i) { return pts[i].y; }
    const double &operator[](int i) const { return pts[i].y; }
};

struct Mesh {
    Region &reg;

    int     nSmooth;
    int     nRegion,nRLine,iRLine; /**/
    int     **boundVert;
    int     *nVert,*region,*bCut;

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
    void regularity(bool Regul);

    XCoord x;
    YCoord y;
public:
    const bool FAF;
    Mesh(Region &reg, const bool FAF) : reg(reg), x(pts), y(pts), FAF(FAF) { }
    Mesh();
    /* exported  function  function */
    int lastpoint() const { return (int)pts.size() - 1; }
    void init( void );
    void addPoint(double x, double y, bool skip_neib = false);
    void addTria( int v1, int v2, int v3, int lab );
    void outMesh( void ) const;
};
#endif
