#ifndef H_TRIA_MESH2D
#define H_TRIA_MESH2D

#include "struct.h"
#include "mesh.h"
#include "tree.h"

#include <set>
#include <vector>

struct edge {
    int v1, v2;
    int &operator[](const int i) {
        if (i == 0)
            return v1;
        if (i == 1)
            return v2;
        throw;
    }
    const int &operator[](const int i) const {
        if (i == 0)
            return v1;
        if (i == 1)
            return v2;
        throw;
    }
    bool operator< (const edge &o) const {
        if (v1 < o.v1)
            return true;
        if (v1 > o.v2)
            return false;
        return v2 < o.v2;
    }
};

class Triangulation {
    double minrho,alpha,beta;
    Mesh &mesh;
    StrucTree2 &tree;

    bool FAF;
    int new_vert;

    PStrucFace2 intedge;
    int ntodelete;
    PStrucFace2 todelete[3];

    std::vector<std::set<int> > eadj;
    std::vector<std::set<edge> > tadj;

    double dist(int  a, int b) const;
    void fill_eadj();
    void fill_tadj();
    int opt_func(int nfixed);

    double func_q(int a, int b, int c) const;
    void func_xy(int c, int a, int b, double dx[2], double delta) const;

    double det2i3(int v1, int v2, int v3) const;
    int idet2i3(int v1, int v2, int v3) const;
    int intsect(int a, int b, int c, int u, int v) const;

    int check(PStrucFace2 e, int pn);
    double height(PStrucFace2 e) const;
    int new_point(PStrucFace2 e);

    int chknadd(int v1, int v2);

    int newTria(int lab);

public:
    Triangulation(Mesh &mesh, StrucTree2 &tree, bool FAF) : mesh(mesh), tree(tree), FAF(FAF) { }
/* error codes:
 *  0 - success
 * -1 - zero sized edge (error in user data)
 * -2 - internal search failed (most likely self intersection in front)
 ****************************************************************************/
    int makeTria();
};

#endif
