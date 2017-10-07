#ifndef H_TRIA_MESH2D
#define H_TRIA_MESH2D

#include "mesh.h"
#include "tree.h"
#include "boundary.h"

#include <set>
#include <vector>

#ifndef SWIG

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

#endif

class Triangulation {
    Mesh mesh;
    Tree tree;
    Boundary bnd;
    const Metric &metric;

    int new_vert;

    PFace intedge;
    int ntodelete;
    PFace todelete[3];

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

    int check(PFace e, int pn);
    double height(PFace e) const;
    int new_point(PFace e);
    int chknadd(int v1, int v2);
    int newTria(int lab);

    double sizeFace(double x, double y) {
        const auto tr = bnd.bb.fromUnit();
        return metric.size(tr(coord(x, y))) / tr.m;
    }

    void meshBoundary();
    void generateInitialFront(int region);

public:
    Triangulation(const Boundary &bnd, const Metric &metric) : bnd(bnd), metric(metric) { }
//    Triangulation(Mesh &mesh, Tree &tree) : mesh(mesh), tree(tree) { }
/* error codes:
 *  0 - success
 * -1 - zero sized edge (error in user data)
 * -2 - internal search failed (most likely self intersection in front)
 ****************************************************************************/
    int generate(bool topological = true, int nSmooth = 5);
    std::vector<coord> points() const {
        std::vector<coord> ret;
        const auto transf = bnd.bb.fromUnit();
        for (const Vertex &v : mesh.pts) {
            ret.push_back(transf(coord(v.x, v.y)));
        }
        return ret;
    }
    const std::vector<Triangle> &triangles() const {
        return mesh.tri;
    }
    std::vector<Edge> edges() const {
        std::vector<Edge> ret;
        for (size_t i = 0; i < bnd.segments.size(); i++) {
            const Segment &s = *bnd.segments[i];
            std::vector<int> seq;
            seq.push_back(s.begCorner);
            std::pair<int, int> meshedSeg = mesh.seg[i];
            for (int j = meshedSeg.first; j < meshedSeg.second; j++)
                seq.push_back(j);
            seq.push_back(s.endCorner);
            for (size_t j = 0; j < seq.size() - 1; j++)
                ret.emplace_back(seq[j], seq[j+1], s.label);
        }
        return ret;
    }
    void saveToVtk(const std::string &fn) {
        mesh.saveVtk(fn);
    }
};

#endif
