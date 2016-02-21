#ifndef H_MESH_MESH2D
#define H_MESH_MESH2D

#include <vector>
#include <string>

struct Vertex {
    double x, y;
    bool remove;
    bool skip_neib;
    bool alreadySwapped;
    std::vector<int> neib;
    std::vector<int> neibTria;
    Vertex(double x, double y, bool skip_neib = false) : x(x), y(y), remove(false), skip_neib(skip_neib), alreadySwapped(false) { }
    Vertex() : remove(false), skip_neib(false), alreadySwapped(false) { }
    void move(double xx, double yy) { x = xx; y = yy; }
};

struct Triangle {
    int v1, v2, v3;
    int label;
    bool remove;
    Triangle(int v1, int v2, int v3, int lab) : v1(v1), v2(v2), v3(v3), label(lab), remove(false) { }
    Triangle() : remove(false) { }
};

struct Edge {
    int v1, v2;
    int label;
    Edge(int v1, int v2, int label)
        : v1(v1), v2(v2), label(label)
    { }
};

struct Mesh {
    std::vector<std::pair<int, int> > seg;
    std::vector<Vertex> pts;
    std::vector<Triangle> tri;

    std::vector<int> vb;
    std::vector<int> ve;
    std::vector<int> tria1;
    std::vector<int> tria2;

    double angle( int v1, int v2, int v ) const;
    void test_quality();
    void smoothing();
    void delVertex(int n);
    void changeVertex(int v1, int v2);
    void delTria(int n);
    void changeTria(int iTria, int v1, int v2, int v3);
    void sortNeigbor(int iVert);
    void calcNeigTria();
    void calcNeigbor();
    void calcEdge();
    void pack();
    void deleteVertex();
    void mergeVertex();
    void splitVertex();
    void swapEdge();
    void optimize(bool topological, int nSmooth);
public:
    Mesh() { }
    /* exported  function  function */
    int lastpoint() const { return (int)pts.size() - 1; }
    void addVertex(double x, double y, bool skip_neib = false);
    void addTria(int v1, int v2, int v3, int lab);
    void saveVtk(const std::string &fn);
};
#endif
