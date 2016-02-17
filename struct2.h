#ifndef H_STRUCT2_MESH2D
#define H_STRUCT2_MESH2D

#include <vector>

struct Point {
    double x, y;
    bool remove;
    bool skip_neib;
    bool alreadySwapped;
    std::vector<int> neib;
    std::vector<int> neibTria;
    Point(double x, double y, bool skip_neib = false) : x(x), y(y), remove(false), skip_neib(skip_neib), alreadySwapped(false) { }
    Point() : remove(false), skip_neib(false), alreadySwapped(false) { }
    void move(double xx, double yy) { x = xx; y = yy; }
};

struct Triangle {
    int v1, v2, v3;
    int label;
    bool remove;
    Triangle(int v1, int v2, int v3, int lab) : v1(v1), v2(v2), v3(v3), label(lab), remove(false) { }
    Triangle() : remove(false) { }
};

typedef struct{
    int     nSmooth;
    int     nRegion,nRLine,iRLine; /**/
    int     **boundVert;
    int     *nVert,*region,*bCut,*nRPoint,*nRTria;

    std::vector<Point> pts;
    std::vector<Triangle> tri;

//    int     **neigTria;
    std::vector<int> vb;
    std::vector<int> ve;
    std::vector<int> tria1;
    std::vector<int> tria2;

    double  size;
    char    outFileName[128];
    char    inFileName[128];
    char    parFileName[128];
    char    debFileName[128];
    char    debug;
}  StrucMesh2;


typedef struct{
    int     v1,v2;     /**/
    int     f;          /* number  of  face       */
    double  x,y,s;    /**/
}  StrucFace2;
typedef  StrucFace2  *PStrucFace2;


typedef struct _MtEC{
    PStrucFace2  face;
    struct _MtEC *next;
} entrychain;

typedef struct _MtSN2{
    int           entrycount;
    entrychain    *firstentry;
    struct _MtSN2 *parent;
    struct _MtSN2 *nodelist[4];
} StrucNode2d;
typedef  StrucNode2d  *PStrucNode2d;

typedef struct {
    int           flag;
    union {
        PStrucFace2 faceNode;
        int v;
    } u;
} StrucNode2;
typedef  StrucNode2  *PStrucNode2;


typedef  PStrucNode2  StrucList2[4];
typedef  StrucList2   *PStrucList2;

#define  MAX_NEIGBOR2   16
#define  S_StrucList2  sizeof(StrucList2)
#define  S_StrucNode2  sizeof(StrucNode2)
#define  S_StrucFace2  sizeof(StrucFace2)

#define  S_StrucNode2d  sizeof(StrucNode2d)
#define  S_entrychain  sizeof(entrychain)


/* exported  function  function */
void init( void );
void addPoint( double x, double y );
void addTria( int v1, int v2, int v3, int lab );
void outMesh( void );


#endif
