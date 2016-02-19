#ifndef H_TREE2_MESH2D
#define H_TREE2_MESH2D

#include "struct.h"

#include <vector>

class Mesh;
class Triangulation;

struct Tree {
    PStrucNode2d root;            /*  root  of  the  quadtree  */

    std::vector<PStrucFace2> faces;
    std::vector<PStrucFace2> vicinityFaces;

    double       xVicinity,yVicinity;   /*  center  of  vicinity */
    double       sVicinity;       /*  size  of  vicinity */

    int          fill,empty;      /*  global  for  recursive  remove  function  */
    double       xc,yc,side;      /*  global  for  recursive  remove  &  insert  function  */
    double       x,y;             /*  global  for  recursive  remove  &  insert  function  */

    void addFaceArray(PStrucFace2 face);
    void remFaceArray(int face);
    int  sectQuad() const;
    void vicinityFacesRec(PStrucNode2d node);
public:
    Tree();
    ~Tree();
    static int direction(double x,double y,double xc,double yc);
    static void center(double *x,double *y,double side,int d);
    static double distance(double x,double y,double xc,double yc);
    static double distanceS(double x,double y,double xc,double yc);
    PStrucFace2 addFace(const Mesh &mesh, int v1, int v2, int twin);
    void remFace(PStrucFace2  face);
    void buildVicinityFaces(double x, double y, double size);

    friend class Triangulation;
};

#endif
