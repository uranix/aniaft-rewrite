#ifndef H_TREE2_MESH2D
#define H_TREE2_MESH2D

#include <vector>

class Mesh;
class Triangulation;

typedef struct Face {
    int     v1,v2;     /**/
    int     f;          /* number  of  face       */
    double  x,y,s;    /**/
} *PFace;

typedef struct _MtEC {
    PFace  face;
    struct _MtEC *next;
} entrychain;

typedef struct _MtSN2{
    int           entrycount;
    entrychain    *firstentry;
    struct _MtSN2 *parent;
    struct _MtSN2 *nodelist[4];
} Node2d;
typedef  Node2d  *PNode2d;

typedef struct {
    int           flag;
    union {
        PFace faceNode;
        int v;
    } u;
} Node2;
typedef  Node2  *PNode2;


typedef  PNode2  List2[4];
typedef  List2   *PList2;

struct Tree {
    PNode2d root;            /*  root  of  the  quadtree  */

    std::vector<PFace> faces;
    std::vector<PFace> vicinityFaces;

    double       xVicinity,yVicinity;   /*  center  of  vicinity */
    double       sVicinity;       /*  size  of  vicinity */

    int          fill,empty;      /*  global  for  recursive  remove  function  */
    double       xc,yc,side;      /*  global  for  recursive  remove  &  insert  function  */
    double       x,y;             /*  global  for  recursive  remove  &  insert  function  */

    void addFaceArray(PFace face);
    void remFaceArray(int face);
    int  sectQuad() const;
    void vicinityFacesRec(PNode2d node);
public:
    Tree();
    ~Tree();
    static int direction(double x,double y,double xc,double yc);
    static void center(double *x,double *y,double side,int d);
    PFace addFace(const Mesh &mesh, int v1, int v2);
    void remFace(PFace  face);
    void buildVicinityFaces(double x, double y, double size);

    friend class Triangulation;
};

#endif
