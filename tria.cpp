#include <cmath>
#include <cstdlib>

#include "region.h"
#include "tria.h"
#include "user.h"

#include <map>

#define SHOWPROGRESS 0

extern Mesh mesh;
extern Tree tree;

double Triangulation::dist(int  a, int b) const {
    const Point &pa = mesh.pts[a];
    const Point &pb = mesh.pts[b];
    return std::sqrt((pa.x - pb.x)*(pa.x - pb.x) + (pa.y - pb.y)*(pa.y - pb.y));
}

double Triangulation::func_q(int a, int b, int c) const {
    const Point &pa = mesh.pts[a];
    const Point &pb = mesh.pts[b];
    const Point &pc = mesh.pts[c];
    const double S = (
            pb.x*pa.y - pa.x*pb.y +
            pc.x*(pb.y-pa.y) +
            pc.y*(pa.x-pb.x)
        )/2.0;
    const double L = 2.0*pc.x*pc.x + 2.0*pa.x*pa.x + 2.0*pb.x*pb.x
        - 2.0*pc.x*pa.x - 2.0*pc.x*pb.x - 2.0*pa.x*pb.x +
        2.0*pc.y*pc.y + 2.0*pa.y*pa.y + 2.0*pb.y*pb.y
        - 2.0*pc.y*pa.y - 2.0*pc.y*pb.y - 2.0*pa.y*pb.y;
    return  4.0*std::sqrt(3.0)*S/L;
}

void Triangulation::func_xy(int c, int a, int b, double dx[2], double delta) const {
    const Point &pa = mesh.pts[a];
    const Point &pb = mesh.pts[b];
    const Point &pc = mesh.pts[c];
    const double S = (
            pb.x*pa.y - pa.x*pb.y +
            pc.x*(pb.y-pa.y) +
            pc.y*(pa.x-pb.x)
        )/2.0;
    const double H = S + std::sqrt(S*S + 4.0*delta*delta);
    //    H = 2.0*S;
    const double L = 2.0*pc.x*pc.x + 2.0*pa.x*pa.x + 2.0*pb.x*pb.x
                   - 2.0*pc.x*pa.x - 2.0*pc.x*pb.x - 2.0*pa.x*pb.x
                   + 2.0*pc.y*pc.y + 2.0*pa.y*pa.y + 2.0*pb.y*pb.y
                   - 2.0*pc.y*pa.y - 2.0*pc.y*pb.y - 2.0*pa.y*pb.y;

    const double Lx = 4.0*pc.x - 2.0*(pa.x+pb.x);
    const double Ly = 4.0*pc.y - 2.0*(pa.y+pb.y);

    const double G = std::sqrt(S*S + 4*delta*delta);
    const double s3 = 1.7320508075688772;

    const double Hx = (pb.y-pa.y)/2.0 + (S*(pb.y-pa.y))/2.0/G;
    const double Hy = (pa.x-pb.x)/2.0 + (S*(pa.x-pb.x))/2.0/G;

    const double q = 2.0*L/H/4.0/s3;
    const double q2 = q * q;
    const double q4 = q2 * q2;
    const double f = q * q2 * q4;

    dx[0] = (2.0*Lx*H - 2.0*Hx*L)/H/H/4.0/s3 * f;
    dx[1] = (2.0*Ly*H - 2.0*Hy*L)/H/H/4.0/s3 * f;
}

static double lldet2(double a, double b, double c, double d) {
    return  a*d - b*c;
}

static double orient2d(double a1, double a2, double b1, double b2, double c1, double c2) {
    double a = a1 - c1,  b = b1 - c1;
    double c = a2 - c2,  d = b2 - c2;
    return lldet2(a, b, c, d);
}

double Triangulation::det2i3(int v1, int v2, int v3) const {
    const Point &p1 = mesh.pts[v1];
    const Point &p2 = mesh.pts[v2];
    const Point &p3 = mesh.pts[v3];
    return orient2d(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
}

int Triangulation::idet2i3(int v1, int v2, int v3) const {
    double d = det2i3(v2, v1, v3);
    if (d > 1e-16) 
        return +1;
    if (d < -1e-16)
        return -1;
    return  0;
}

int Triangulation::intsect(int a, int b, int c, int u, int v) const {
    int dup = 0;

    if ((u == a) || (u == b) || (u == c))
        dup++;
    if ((v == a) || (v == b) || (v == c))
        dup++;
    if (dup == 0) {
        int uv = idet2i3(u, v, a) + idet2i3(u, v, b) + idet2i3(u, v, c);
        if ((uv == 3) || (uv == -3)) return 0;
        if (idet2i3(b, c, u) + idet2i3(b, c, v) == -2) return  0;
        if (idet2i3(c, a, u) + idet2i3(c, a, v) == -2) return  0;
        if (idet2i3(a, b, u) + idet2i3(a, b, v) == -2) return  0;
    } else if (dup == 1) {
        if (idet2i3(b, c, u) + idet2i3(b, c, v) == -1) return  0;
        if (idet2i3(c, a, u) + idet2i3(c, a, v) == -1) return  0;
        if (idet2i3(a, b, u) + idet2i3(a, b, v) == -1) return  0;
    } else
        return  0;
    return  1;
}

int Triangulation::check(PStrucFace2 e, int pn) {
    int v1, v2, p1, p2;

    intedge = NULL;
    v1 = e->v1;
    v2 = e->v2;

    if (idet2i3(v1, v2, pn) != 1) {
        //fprintf(stderr, "\ninverted? %d (%d %d %d)\n", idet2i3(v1, v2, pn), v1, v2, pn);
        return 1;
    }

    for (const PStrucFace2 face : tree.vicinityFaces) {
        p1 = face->v1;
        p2 = face->v2;
        if (intsect(v1, v2, pn, p1, p2)) {
            intedge = face;
            return 1;
        }
    }
    return 0;
}

double Triangulation::height(PStrucFace2 e) const {
    int v0 = e->v1;
    int v1 = e->v2;

    const Point &np = mesh.pts[new_vert];
    const Point &p0 = mesh.pts[v0];
    const Point &p1 = mesh.pts[v1];
    double x  = np.x;
    double y  = np.y;
    double x0 = p0.x;
    double y0 = p0.y;
    double x1 = p1.x;
    double y1 = p1.y;

    double c = ((x-x0)*(x1-x0) + (y-y0)*(y1-y0)) / ((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
    if (c<=0.0)
        return std::sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
    else if (c>=1.0)
        return std::sqrt((x1-x)*(x1-x) + (y1-y)*(y1-y));

    double s = (x0-x)*(y1-y) - (y0-y)*(x1-x);
    double r = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));

    return std::fabs(s) / r;
}

#define CND_MAX 1000

struct BestCandidates {
    size_t cnd_max;

    std::vector<std::pair<int, double> > cnd;

    BestCandidates(size_t cnd_max = CND_MAX) : cnd_max(cnd_max) {
    }

    bool insert(int v, double q) {
        decltype(cnd)::iterator k;
        for (k = cnd.begin(); k != cnd.end(); k++) {
            if (k->first == v)
                return false;
            if (q > k->second)
                break;
        }
        cnd.insert(k, std::make_pair(v, q));
        if (cnd.size() > cnd_max)
            cnd.pop_back();
        return true;
    }
};

// called from newTria
int Triangulation::new_point(PStrucFace2 e) {
    int nchk, chk[CND_MAX];

    int v1 = e->v1;
    int v2 = e->v2;

    double x1 = mesh.pts[v1].x;
    double y1 = mesh.pts[v1].y;
    double x2 = mesh.pts[v2].x;
    double y2 = mesh.pts[v2].y;

    double b = x1 - x2;
    double a = y2 - y1;
    double p = std::sqrt(a * a + b * b);

    if (p == 0.0)
        return -1;
    a /= p;
    b /= p;
    double r = p / 2.0;

    double x, y;

    if (!mesh.FAF) {
        x = 0.5 * (x1 + x2) + a * 0.3 * p;
        y = 0.5 * (y1 + y2) + b * 0.3 * p;
        p = sizeFace(x, y);
        if (p*p-r*r < p*p*3.0/4.0)
            p = std::sqrt(r*r + p*p*3.0/4.0);
    } else {
        const double COARSEFACTOR = 1.5;
        p = COARSEFACTOR * std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    }
    if (p < 1.1*r)
        p = 1.1*r;

    x = 0.5 * (x1 + x2) + a * std::sqrt(p*p-r*r);
    y = 0.5 * (y1 + y2) + b * std::sqrt(p*p-r*r);

    mesh.addPoint(x, y);
    new_vert = mesh.pts.size() - 1;

    /*____________________ TEST __________ TEST ____________________________*/

    double radius = Tree::distance(x, y, x1, y1);

    r = radius * 1.0001220703125; // 8193. / 8192
    double rmin = (beta*r + minrho + std::sqrt((beta*r - minrho)*(beta*r - minrho) + alpha))/2.0;
    r = radius * 2.0;

    tree.buildVicinityFaces(x, y, r);

    int neari = -1;
    double rv = rmin;
    double hc = height(e);
    int dirty = 0;

    BestCandidates cand;
    cand.insert(new_vert, func_q(v1, v2, new_vert));

    for (const PStrucFace2 face : tree.vicinityFaces) {
        double h = height(face);
        if (h < 0.5*hc) {
            dirty++;
//            printf("dirty: %lf, %lf\n", hc, h);
        }
        for (const int pn : { face->v1, face->v2 } ) {
            if ((pn==v1) || (pn==v2))
                continue;
            if (idet2i3(v1, v2, pn) != 1)
                continue;
            r = Tree::distance(mesh.pts[pn].x, mesh.pts[pn].y, x, y);
            if (r < rv) {
                rv = r;
                neari = pn;
            }
            double q = func_q(v1, v2, pn);

            cand.insert(pn, q);
        }
    }

    nchk = 0;
    if (dirty && cand.cnd.size() > 1)  {
        chk[nchk++] = new_vert;
        new_vert = cand.cnd[1].first;
    }
    if (neari >= 0)
        new_vert = neari;

    for ( ; nchk < CND_MAX; nchk++) {
        chk[nchk] = new_vert;
        if (check(e, new_vert) == 0)
            return 0;

        bool broke = false;
        for (size_t k = 0; k < cand.cnd.size(); k++) {
            new_vert = cand.cnd[k].first;
            // search for new_vert in chk
            int m;
            for (m=0; m <= nchk; m++) {
                if (new_vert == chk[m]) {
                    broke = true;
                    break;
                }
            }
            if (m>nchk) {
                broke = true;
                break;
            }
        }
        if (!broke) {
            return -2;
        }

    }
    return  1;
}

int Triangulation::chknadd(int v1, int v2) {
    int p1, p2;

    for (const PStrucFace2 face : tree.vicinityFaces) {
        p1 = face->v1;
        p2 = face->v2;
        if ((p1==v2) && (p2==v1)) {
            todelete[ntodelete++] = face;
            return 1;
        }
    }
    tree.addFace(mesh, v1, v2, 0);
    return 0;
}

int Triangulation::newTria(int lab) {
    int nn;
    int v1, v2;
    PStrucFace2 e1;

    e1 = tree.faces.front();
    v1 = e1->v1;
    v2 = e1->v2;
    nn = new_point(e1);
    if (nn != 0)
        return nn;

    if (new_vert != (int)mesh.pts.size() - 1)
        mesh.pts.pop_back();

    todelete[0] = e1;
    ntodelete = 1;
    mesh.addTria(v1, v2, new_vert, lab);

    chknadd(new_vert, v2);
    chknadd(v1, new_vert);
    for (nn = 0; nn < ntodelete; nn++)
        tree.remFace(todelete[nn]);
    return 0;
}

extern int    nVRTglobal;
extern int    nTRIglobal;

int Triangulation::makeTria() {
    int i = 0, err = 0, smooth = 0;

    initRegion();

    for (i = 1; i <= mesh.nRegion; i++) {
        int nRPointPrev = mesh.pts.size();
        minrho = dist(tree.faces[0]->v1, tree.faces[0]->v2);

        for (const auto &f : tree.faces)
            if (minrho > dist(f->v1, f->v2))
                minrho = dist(f->v1, f->v2);

        if (mesh.FAF)
            minrho *= 0.95;
        else
            minrho = 0.0;
        beta = 0.5;
        minrho *= beta;
        alpha = 2.0*minrho*(1.0-beta)/beta;
        alpha *= alpha;

        while (!tree.faces.empty()) {
            err = newTria(i);
            if (err)
                break;
            if (SHOWPROGRESS && (mesh.tri.size() % 100 == 0)) {
                printf("Number of Point = %lu Number of Tria = %lu\n", mesh.pts.size(), mesh.tri.size());
                fflush(stdout);
            }
        }

        if (err)
            break;
        if (SHOWPROGRESS)
            printf("\n");
        fill_eadj();
        fill_tadj();

        smooth += opt_func(nRPointPrev);

        if (i != mesh.nRegion) {
            initAddRegion(i + 1);
        }
    }

    printf("\nRESULT :  %5lu     %5lu    \n", mesh.pts.size(), mesh.tri.size());

    if (!smooth)
        mesh.regularity(true);

    mesh.test_quality();
    mesh.outMesh();

    if (!err) {
        nVRTglobal = mesh.pts.size();
        nTRIglobal = mesh.tri.size();
    }

    return err;
}

void Triangulation::fill_eadj() {
    eadj.clear();
    eadj.resize(mesh.pts.size());

    for (const auto t : mesh.tri) {
        int a, b, c;
        a = t.v1;
        b = t.v2;
        c = t.v3;

        eadj[a].insert(b);
        eadj[a].insert(c);
        eadj[b].insert(c);
        eadj[b].insert(a);
        eadj[c].insert(a);
        eadj[c].insert(b);
    }
}

void Triangulation::fill_tadj() {
    tadj.clear();
    tadj.resize(mesh.pts.size());

    for (const auto t : mesh.tri) {
        int a;
        edge b;

        a    = t.v1;
        b[0] = t.v2;
        b[1] = t.v3;

        tadj[a].insert(b);

        a    = t.v2;
        b[0] = t.v3;
        b[1] = t.v1;
        tadj[a].insert(b);

        a    = t.v3;
        b[0] = t.v1;
        b[1] = t.v2;
        tadj[a].insert(b);
    }
}

int Triangulation::opt_func(int nfixed) {
    double x[2], z[2];

    std::vector<int> ps;

    double ds = 0.0;
    for (size_t p = nfixed; p < mesh.pts.size(); p++) {
        ps.push_back(p);

        double d = 0.0;
        int s = 0;
        for (const int &c : eadj[p]) {
            d += dist(p, c);
            s++;
        }
        if (s > 0)
            d /= s;
        ds += d;
    }
    if (ps.size() == 0)
        return 0;

    ds /= ps.size();
    //printf("ds=%lf\n", ds);

    typedef double vertex[2];
    std::vector<vertex> dx(ps.size());
    std::vector<Point> bkp(ps.size());

    for (size_t j = 0; j < ps.size(); j++) {
        int p = ps[j];
        bkp[j] = mesh.pts[p];
    }

    const int n_iters = 100; // 4000
    double d = 0.25 * ds * ds / n_iters;  // 1.0

    for (int s = 0; s < n_iters; s++) {
        const double delta = 0.01 * ds * ds * (1.0 - 0.9*s/n_iters);

        //	delta = -as / n_iters;
        double rs = 0.0;
        for (size_t j = 0; j < ps.size(); j++) {
            int p = ps[j];
            x[0] = 0.0;
            x[1] = 0.0;
            for (const edge &e : tadj[p]) {
                int b = e[0];
                int c = e[1];
                func_xy(p, b, c, z, delta);
                x[0] -= z[0];
                x[1] -= z[1];
            }
            double r = std::sqrt(x[0]*x[0] + x[1]*x[1]);
            if (r > 1.0/ds) {
                r = 1.0/ds*(1.0 + log(r*ds))/r;
                x[0] *= r;
                x[1] *= r;
            }
            dx[j][0] = x[0];
            dx[j][1] = x[1];
            r = std::sqrt(x[0]*x[0] + x[1]*x[1]);
            if (rs < r)
                rs = r;
        }
        for (size_t j = 0; j < ps.size(); j++) {
            int p = ps[j];
            mesh.pts[p].x += d*dx[j][0];
            mesh.pts[p].y += d*dx[j][1];
        }
        if (SHOWPROGRESS) {
            printf(" Smoothing %d%%\n", s+1);
            fflush(stdout);
        }
    }
    if (SHOWPROGRESS)
        printf(" Smoothing done\n");

    int s = 0;
    for (const auto &t : mesh.tri)
        if (func_q(t.v1, t.v2, t.v3) <= 0.0)
            s++;

    if (s) {
        fprintf(stderr,"aniAFT: Quality improvement failed, falling back to simple smoothing\n");
        for (size_t j = 0; j < ps.size(); j++) {
            int p = ps[j];
            mesh.pts[p] = bkp[j];
        }
    }

    return s;
}
