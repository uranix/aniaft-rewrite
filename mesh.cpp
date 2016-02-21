#include "mesh.h"
#include "tree.h"
#include "metric.h"

#include <cmath>
#include <cassert>
#include <map>
#include <algorithm>
#include <string>
#include <fstream>

void Mesh::addVertex(double x, double y, bool skip_neib) {
    pts.push_back(Vertex(x, y, skip_neib));
}

void Mesh::addTria(int v1, int v2, int v3, int lab) {
    tri.push_back(Triangle(v1, v2, v3, lab));
}

void Mesh::smoothing() {
    for (Vertex &p : pts) {
        double xx=0.;
        double yy=0.;
        for (const auto &nei : p.neib) {
            double x = pts[nei].x;
            double y = pts[nei].y;
            xx += x;
            yy += y;
        }
        int n = p.neib.size();
        if( n > 0 ){
            xx /= n;
            yy /= n;
            p.move(xx, yy);
        }
    }
}

void Mesh::delVertex(int n) {
    pts[n].remove = true;
}

void Mesh::changeVertex(int v1, int v2) {
    pts[v1].x = 0.5*(pts[v1].x + pts[v2].x);
    pts[v1].y = 0.5*(pts[v1].y + pts[v2].y);
}

void Mesh::delTria(int n) {
    tri[n].remove = true;
}

void Mesh::changeTria(int iTria, int v1, int v2, int v3) {
    tri[iTria].v1 = v1;
    tri[iTria].v2 = v2;
    tri[iTria].v3 = v3;
}

void Mesh::sortNeigbor(int iVert) {
    std::vector<std::pair<int, double> > na;

    double x0 = pts[iVert].x;
    double y0 = pts[iVert].y;

    for (const int nei : pts[iVert].neib) {
        double x = pts[nei].x - x0;
        double y = pts[nei].y - y0;
        double p = sqrt(x*x+y*y);
        if (p == 0.)
            p = 1;
        x /= p;
        y /= p;
        double angle;

        if (x < 0.)
            angle = 2. - y;
        else if (y >= 0.)
            angle = y;
        else
            angle = 4. + y;
        na.push_back(std::make_pair(nei, angle));
    }

    std::sort(na.begin(), na.end(), []
        (const std::pair<int, double> &a, const std::pair<int, double> &b) {
            return a.second < b.second;
        }
    );

    std::transform(na.begin(), na.end(), pts[iVert].neib.begin(),
            [] (const std::pair<int, double> &v) { return v.first; });
}

void Mesh::calcNeigTria() {
    for (Vertex &p : pts) {
        p.neibTria.clear();
        p.alreadySwapped = false;
    }
    int i = 0;
    for (const auto tr : tri) {
        pts[tr.v1].neibTria.push_back(i);
        pts[tr.v2].neibTria.push_back(i);
        pts[tr.v3].neibTria.push_back(i);
        i++;
    }
}

void Mesh::calcNeigbor() {
    std::vector<int> vert;

    for (Vertex &p : pts) {
        if (!p.skip_neib)
            p.neib.clear();
    }

    for (int i = 0; i < (int)pts.size(); i++) {
        Vertex &p = pts[i];
        if (p.skip_neib)
            continue;
        vert.clear();
        for (const int itr : p.neibTria) {
            const Triangle &tr = tri[itr];
            vert.push_back(tr.v1);
            vert.push_back(tr.v2);
            vert.push_back(tr.v3);
        }
        for (const int v : vert) {
            if (v == i)
                continue;

            bool already = false;
            for (const int &nv : p.neib)
                if (v == nv)
                    already = true;
            if (already)
                continue;
            p.neib.push_back(v);
        }
    }
}

void Mesh::calcEdge() {
    calcNeigTria();
    calcNeigbor();

    vb.clear();
    ve.clear();
    tria1.clear();
    tria2.clear();

    for(int i = 0; i < (int)pts.size(); i++) {
        Vertex &p = pts[i];
        if (p.skip_neib)
            continue;
        for (const int vert : p.neib) {
            if (vert <= i)
                continue;
            if (pts[vert].skip_neib)
                continue;

            vb.push_back(i);
            ve.push_back(vert);
            bool swap = true;
            for (size_t itr : p.neibTria) {
                const Triangle &tr = tri[itr];
                if (tr.v1 == vert || tr.v2 == vert || tr.v3 == vert) {
                    if (swap) {
                        tria1.push_back(itr);
                        swap = false;
                    } else
                        tria2.push_back(itr);
                }
            }
            swap = true;
            int itr = tria1.back();
            const Triangle &tr = tri[itr];
            if (tr.v1 == i && tr.v2 == vert)
                swap = false;
            else if (tr.v2 == i && tr.v3 == vert)
                swap = false;
            else if (tr.v3 == i && tr.v1 == vert)
                swap = false;
            if (swap)
                std::swap(tria1.back(), tria2.back());
        }
    }
}

template<class Cont, class Pred>
void drop(Cont &cont, Pred pred) {
    cont.erase(std::remove_if(cont.begin(), cont.end(), pred), cont.end());
}

void Mesh::pack() {
    drop(tri, [] (const Triangle &t) { return t.remove; });

    std::vector<int> newidx;

    int idx = 0;
    for (const Vertex &p : pts) {
        if (p.remove)
            newidx.push_back(-1);
        else {
            newidx.push_back(idx);
            idx++;
        }
    }

    for (auto &tr : tri) {
        tr.v1 = newidx[tr.v1];
        tr.v2 = newidx[tr.v2];
        tr.v3 = newidx[tr.v3];
        assert(tr.v1 >= 0 && tr.v2 >= 0 && tr.v3 >= 0);
    }

    drop(pts, [] (const Vertex &p) { return p.remove; });
}

// Deletes points with 3 or 4 neighbors
void Mesh::deleteVertex() {
    for(int i = 0; i < (int)pts.size(); i++) {
        Vertex &p = pts[i];
        int j = p.neib.size();
        if (j == 3) {
            delVertex(i);
            sortNeigbor(i);
            changeTria(p.neibTria[0], p.neib[1], p.neib[0], p.neib[2]);
            delTria   (p.neibTria[1]);
            delTria   (p.neibTria[2]);
        }
        else if( j == 4 ){
            delVertex(i);
            sortNeigbor(i);
            int n1 = pts[p.neib[0]].neib.size() - 7;
            int n2 = pts[p.neib[1]].neib.size() - 7;
            int n3 = pts[p.neib[2]].neib.size() - 7;
            int n4 = pts[p.neib[3]].neib.size() - 7;
            int sum  = (n1+1)*(n1+1) + (n3+1)*(n3+1) + n2*n2 + n4*n4;
            int swap = (n2+1)*(n2+1) + (n4+1)*(n4+1) + n1*n1 + n3*n3;
            if (swap < sum) {
                changeTria(p.neibTria[0], p.neib[1], p.neib[0], p.neib[3]);
                changeTria(p.neibTria[1], p.neib[3], p.neib[2], p.neib[1]);
            } else {
                changeTria(p.neibTria[0], p.neib[1], p.neib[0], p.neib[2]);
                changeTria(p.neibTria[1], p.neib[3], p.neib[2], p.neib[0]);
            }
            delTria(p.neibTria[2]);
            delTria(p.neibTria[3]);
        }
    }
    pack();
    calcNeigTria();
    calcNeigbor();
}

// Merges two points with 5 neibs
void Mesh::mergeVertex() {
    calcEdge();
    for(size_t i = 0; i < vb.size(); i++){
        int v1 = vb[i];
        int v2 = ve[i];
        int n1 = pts[v1].neib.size();
        int n2 = pts[v2].neib.size();
        if( n1 == 5 && n2 == 5 ){
            // Dont merge futher these points
            pts[v1].neib.push_back(-1);
            pts[v2].neib.push_back(-1);
            changeVertex(v1, v2);
            delVertex(v2);
            for (const int itr : pts[v2].neibTria) {
                Triangle &tr = tri[itr];
                if( tr.v1 == v2 )
                    tr.v1 = v1;
                if( tr.v2 == v2 )
                    tr.v2 = v1;
                if( tr.v3 == v2 )
                    tr.v3 = v1;
            }
            delTria(tria1[i]);
            delTria(tria2[i]);
        }
    }

    pack();
    calcNeigTria();
    calcNeigbor();
}

// splits points with more than 7 neibs
void Mesh::splitVertex() {
    for (int i = 0; i < (int)pts.size(); i++){
        Vertex &p = pts[i];
        int n1 = p.neib.size();
        if( n1 <= 7 )
            continue;
        int n2 = n1/2;
        n1 = n1 - n2;
        sortNeigbor(i);

        const Vertex &p2 = pts[p.neib[0]];

//        double size = sizeFace(p.x, p.y);
        double size = Metric::distance(p.x, p.y, p2.x, p2.y);

        addVertex(p.x+0.1*size, p.y);
        Vertex &pv = pts.back();
        int vert = lastpoint();

        pv.neib.push_back(-1);

        for (const int nei : p.neib) {
            Vertex &pn = pts[nei];
            if (!pn.skip_neib) {
                pn.neib.clear();
                pn.neib.push_back(-1);
            }
        }

        int lab = -1;
        for(int j = 0; j < n1 - 1; j++) {
            Triangle &tr = tri[p.neibTria[j]];
            if (lab == -1)
                lab = tr.label;

            if (lab != tr.label)
                throw std::logic_error("aniAFT: different materials inside region");
            changeTria(p.neibTria[j], p.neib[j+1], p.neib[j], i);
        }
        for(int j = n1; j < n1 + n2 - 1; j++) {
            Triangle &tr = tri[p.neibTria[j]];
            if (lab != tr.label)
                throw std::logic_error("aniAFT: different materials inside region");
            changeTria(p.neibTria[j], p.neib[j+1], p.neib[j], vert);
        }
        {
            int it = p.neibTria[n1 - 1];
            Triangle &tr = tri[it];
            if (lab != tr.label)
                throw std::logic_error("aniAFT: different materials inside region");

            changeTria(it, p.neib[n1-1], i, vert);
        }
        {
            int it = p.neibTria[n1 + n2 - 1];
            Triangle &tr = tri[it];
            if (lab != tr.label)
                throw std::logic_error("aniAFT: different materials inside region");
            changeTria(it, p.neib[n1], p.neib[n1-1], vert);
        }
        addTria(p.neib[n1+n2-1], vert, i, lab);
        addTria(p.neib[n1+n2-1], i, p.neib[0],lab);
    }
    calcNeigTria();
    calcNeigbor();

    for(int i=0; i<5; i++)
        smoothing();
}

void Mesh::swapEdge() {
    calcEdge();
    for(size_t i = 0; i < vb.size(); i++) {
        int vb = this->vb[i];
        int ve = this->ve[i];
        int tria1 = this->tria1[i];
        int tria2 = this->tria2[i];

        const Triangle &tr1 = tri[tria1];
        const Triangle &tr2 = tri[tria2];

        // v1 = tr1.v1 ^ tr1.v2 ^ tr1.v3 ^ vb ^ ve
        int v1 = tr1.v1;
        if (v1 == vb || v1 == ve) {
            v1 = tr1.v2;
            if (v1 == vb || v1 == ve)
                v1 = tr1.v3;
        }

        // v2 = tr2.v1 ^ tr2.v2 ^ tr2.v3 ^ vb ^ ve
        int v2 = tr2.v1;
        if (v2 == vb || v2 == ve) {
            v2 = tr2.v2;
            if (v2 == vb || v2 == ve)
                v2 = tr2.v3;
        }

        if (pts[v1].skip_neib || pts[v2].skip_neib)
            continue;
        if (pts[vb].alreadySwapped && pts[ve].alreadySwapped)
            continue;
        int nb = pts[vb].neib.size() - 6;
        int ne = pts[ve].neib.size() - 6;
        int n1 = pts[v1].neib.size() - 6;
        int n2 = pts[v2].neib.size() - 6;
        if (nb < 1 && ne < 1 && n1 < 1 && n2 < 1)
            continue;
        int sum = nb*nb + ne*ne + n1*n1 + n2*n2;
        int swap = (nb-1)*(nb-1) + (ne-1)*(ne-1) + (n1+1)*(n1+1) + (n2+1)*(n2+1);

        if (swap < sum) {
            pts[vb].neib.pop_back();
            pts[ve].neib.pop_back();
            pts[v1].neib.push_back(-1);
            pts[v2].neib.push_back(-1);
            pts[vb].alreadySwapped = true;
            pts[ve].alreadySwapped = true;
            pts[v1].alreadySwapped = true;
            pts[v2].alreadySwapped = true;
            changeTria(tria1, v2, v1, vb);
            changeTria(tria2, v1, v2, ve);
        }
    }

    calcNeigTria();
    calcNeigbor();
}

void Mesh::optimize(bool topological, int nSmooth) {
    int i;

    if (topological) {
        calcNeigTria();
        calcNeigbor();

        for (int i = 0; i < 5; i++)
            smoothing();

        deleteVertex();
        for(i=0;i<5;i++)
            smoothing();

        mergeVertex();
        mergeVertex();
        mergeVertex();
        for(i=0;i<5;i++)
            smoothing();

        splitVertex();
        splitVertex();
        splitVertex();

        swapEdge();
        for(i=0;i<5;i++)
            smoothing();
        swapEdge();
        for(i=0;i<5;i++)
            smoothing();
        swapEdge();
        for(i=0;i<5;i++)
            smoothing();
        swapEdge();
    }

    for(i=0;i<nSmooth;i++)
        smoothing();
}

void Mesh::saveVtk(const std::string &fn) {
    std::ofstream f(fn, std::ios::binary);
    if (!f)
        return;

    f << "# vtk DataFile Version 3.0\n";
    f << "Mesh\n";
    f << "ASCII\n";
    f << "DATASET UNSTRUCTURED_GRID\n";
    f << "POINTS " << pts.size() << " float\n";
    for (const Vertex &v : pts)
        f << v.x << " " << v.y << " 0\n";
    f << "CELLS " << tri.size() << " " << 4 * tri.size() << "\n";
    for (const Triangle &t : tri)
        f << "3 " << t.v1 << " " << t.v2 << " " << t.v3 << "\n";
    f << "CELL_TYPES " << tri.size() << "\n";
    for (size_t k = 0; k < tri.size(); k++)
        f << "5\n";
    f << "CELL_DATA " << tri.size() << "\n";
    f << "SCALARS region int\nLOOKUP_TABLE default\n";
    for (const Triangle &t : tri)
        f << t.label << "\n";
}
