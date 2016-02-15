#include <cmath>
#include <cstdlib>

#include"region2.h"
#include"refine2.h"
#include"tree2.h"
#include"tria2.h"
#include"user2.h"

#include <set>
#include <map>
#include <vector>

#define SHOWPROGRESS 1

extern  StrucMesh2  mesh2;
extern  StrucTree2  tree2;

extern  int boolFAF;
static int new_vert;

static double minrho,alpha,beta;

static PStrucFace2 intedge;
static int ntodelete;
static PStrucFace2 todelete[3];

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

typedef double vertex[3];

std::vector<std::set<int> > eadj;
std::vector<std::set<edge> > tadj;

static double dist(int  a, int b) {
    const Point &pa = mesh2.pts[a];
    const Point &pb = mesh2.pts[b];
    return std::sqrt((pa.x - pb.x)*(pa.x - pb.x) + (pa.y - pb.y)*(pa.y - pb.y));
}

static double func_q(int a, int b, int c) {
    const Point &pa = mesh2.pts[a];
    const Point &pb = mesh2.pts[b];
    const Point &pc = mesh2.pts[c];
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

static int func_xy(int c, int a, int b, double dx[2], double delta, int n) {
    const Point &pa = mesh2.pts[a];
    const Point &pb = mesh2.pts[b];
    const Point &pc = mesh2.pts[c];
    const double S = (
            pb.x*pa.y - pa.x*pb.y +
            pc.x*(pb.y-pa.y) +
            pc.y*(pa.x-pb.x)
        )/2.0;
    const double H = S + std::sqrt(S*S + 4.0*delta*delta);
    //    H = 2.0*S;
    const double L = 2.0*pc.x*pc.x + 2.0*pa.x*pa.x + 2.0*pb.x*pb.x
        - 2.0*pc.x*pa.x - 2.0*pc.x*pb.x - 2.0*pa.x*pb.x +
        2.0*pc.y*pc.y + 2.0*pa.y*pa.y + 2.0*pb.y*pb.y
        - 2.0*pc.y*pa.y - 2.0*pc.y*pb.y - 2.0*pa.y*pb.y;

    const double Lx = 4.0*pc.x - 2.0*(pa.x+pb.x);
    const double Ly = 4.0*pc.y - 2.0*(pa.y+pb.y);

    const double Hx = (pb.y-pa.y)/2.0 + (S*(pb.y-pa.y))/2.0/std::sqrt(S*S + 4.0*delta*delta);
    const double Hy = (pa.x-pb.x)/2.0 + (S*(pa.x-pb.x))/2.0/std::sqrt(S*S + 4.0*delta*delta);

    const double f = std::pow(2.0*L/H/4.0/std::sqrt(3.0), n - 1);
    dx[0] = (2.0*Lx*H - 2.0*Hx*L)/H/H/4.0/std::sqrt(3.0);
    dx[1] = (2.0*Ly*H - 2.0*Hy*L)/H/H/4.0/std::sqrt(3.0);

    dx[0] *= f;
    dx[1] *= f;

    return 0;
}

static double lldet2(double a, double b, double c, double d) {
    return  a*d - b*c;
}
static double orient2d(double a1, double a2, double b1, double b2, double c1, double c2) {
    double a = a1 - c1,  b = b1 - c1;
    double c = a2 - c2,  d = b2 - c2;
    return lldet2(a, b, c, d);
}
static double det2i3(int v1, int v2, int v3) {
    double r;
    const Point &p1 = mesh2.pts[v1];
    const Point &p2 = mesh2.pts[v2];
    const Point &p3 = mesh2.pts[v3];
    r = orient2d(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
    return  r;
}
static int idet2i3(int v1, int v2, int v3) {
    double d = det2i3(v2, v1, v3);
    if (d > 1e-16) return  +1;
    else if (d < -1e-16) return  -1;
    else return  0;
}

static int intsect(int a, int b, int c, int u, int v) {
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
static int check(PStrucFace2 e, int pn) {
    int v1, v2, i, p1, p2;

    intedge = NULL;
    v1 = e->v1;
    v2 = e->v2;

    if (idet2i3(v1, v2, pn) != 1) {
        //fprintf(stderr, "\ninverted? %d (%d %d %d)\n", idet2i3(v1, v2, pn), v1, v2, pn);
        return 1;
    }

    for (i = 0; i < tree2.nVicinityFace; i++) {
        p1 = tree2.vicinityFace[i]->v1;
        p2 = tree2.vicinityFace[i]->v2;
        if (intsect(v1, v2, pn, p1, p2)) {
            intedge = tree2.vicinityFace[i];
            return  1;
        }
    }
    return  0;
}

static double height(PStrucFace2 e) {
    int v0, v1;
    double x, y, x0, y0, x1, y1, s, r, c;

    v0 = e->v1;
    v1 = e->v2;

    const Point &np = mesh2.pts[new_vert];
    const Point &p0 = mesh2.pts[v0];
    const Point &p1 = mesh2.pts[v1];
    x  = np.x;
    y  = np.y;
    x0 = p0.x;
    y0 = p0.y;
    x1 = p1.x;
    y1 = p1.y;

    c = ((x-x0)*(x1-x0) + (y-y0)*(y1-y0)) / ((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
    if (c<=0.0)
        return std::sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
    else if (c>=1.0)
        return std::sqrt((x1-x)*(x1-x) + (y1-y)*(y1-y));

    s = (x0-x)*(y1-y) - (y0-y)*(x1-x);
    r = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));

    return std::fabs(s)/r;
}

#define CND_MAX 1000
static int new_point(PStrucFace2 e) {
    double p, a, b, r, q;
    double x, y, x1, y1, x2, y2;
    int i, j, k, m, pn;
    int v1, v2;
    int ncnd, cnd[CND_MAX];
    double cnd_q[CND_MAX];
    int nchk, chk[CND_MAX];
    int neari;
    double radius, rmin, rv;
    double hc, h;
    int dirty;

    v1 = e->v1;
    v2 = e->v2;

    x1 = mesh2.pts[v1].x;
    y1 = mesh2.pts[v1].y;
    x2 = mesh2.pts[v2].x;
    y2 = mesh2.pts[v2].y;

    b = x1 - x2;
    a = y2 - y1;
    p = std::sqrt(a * a + b * b);

    if (p == 0.0)
        return -1;
    a /= p;
    b /= p;
    r = p / 2.0;

    if (!boolFAF) {
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

    addPoint(x, y);
    mesh2.nPoint--;
    new_vert = mesh2.nPoint;

    /*____________________ TEST __________ TEST ____________________________*/

    radius = distance(x, y, x1, y1);

    r = radius * 1.0001220703125; // 8193. / 8192
    rmin = (beta*r + minrho + std::sqrt((beta*r - minrho)*(beta*r - minrho) + alpha))/2.0;
    r = radius * 2.0;

    vicinityFaces(x, y, r);

    struct {
        int v;
        double qual;
    } candidate;

    std::vector<candidate> canidates;

    candidates.push_back(candidate{new_vert, func_q(v1, v2, new_vert)});
/*
    cnd[0] = new_vert;
    ncnd = 1;
    cnd_q[0] = func_q(v1, v2, new_vert);
    */

    nchk = 0;
    neari = -1;
    rv = rmin;
    hc = height(e);
    dirty = 0;
    for (i = 0; i < tree2.nVicinityFace; i++) {
        h = height(tree2.vicinityFace[i]);
        if (h < 0.5*hc)
            dirty++;//,  printf("dirty: %lf, %lf\n", hc, h);
        for (j = 0, pn = tree2.vicinityFace[i]->v1; j < 2; j++, pn = tree2.vicinityFace[i]->v2) {
            if ((pn==v1) || (pn==v2))
                continue;
            if (idet2i3(v1, v2, pn) != 1)
                continue;
            r = distance(mesh2.x[pn], mesh2.y[pn], x, y);
            if (r < rv) {
                rv = r;
                neari = pn;
            }
            q = func_q(v1, v2, pn);
            for (k=0; k<ncnd; k++) {
                if (pn == cnd[k])
                    break;
                if (q > cnd_q[k])
                    break;
            }
            if (k>=CND_MAX)
                continue;
            if (k==ncnd) {
                cnd[ncnd] = pn;
                cnd_q[ncnd] = q;
                ncnd++;
            } else {
                if (pn == cnd[k])
                    continue;
                for (m=ncnd; m>k; m--) {
                    cnd[m] = cnd[m-1];
                    cnd_q[m] = cnd_q[m-1];
                }
                cnd[k] = pn;
                cnd_q[k] = q;
                ncnd++;
            }
        }
    }
    if (dirty && ncnd > 1)  {
        chk[nchk++] = new_vert;
        new_vert = cnd[1];
    }
    if (neari >= 0)
        new_vert = neari;

    for ( ; nchk < CND_MAX; nchk++) {
        chk[nchk] = new_vert;
        if (check(e, new_vert) == 0)
            return 0;
        for (k=0; k < ncnd; k++) {
            new_vert = cnd[k];
            for (m=0; m <= nchk; m++) {
                if (new_vert == chk[m])
                    break;
            }
            if (m>nchk)
                break;
        }
        if (k>=ncnd) {
            return -2;
        }
    }
    return  1;
} /* new_point */

static int chknadd(int v1, int v2) {
    int i, p1, p2;

    for (i = 0; i < tree2.nVicinityFace; i++) {
        p1 = tree2.vicinityFace[i]->v1;
        p2 = tree2.vicinityFace[i]->v2;
        if ((p1==v2) && (p2==v1)) {
            todelete[ntodelete++] = tree2.vicinityFace[i];
            return 1;
        }
    }
    addFace(v1, v2, 0);
    return 0;
}

static int newTria(int lab) {
    int nn;
    int v1, v2;
    PStrucFace2 e1;

    e1 = tree2.face[0];
    v1 = e1->v1;
    v2 = e1->v2;
    nn = new_point(e1);
    if (nn != 0)
        return nn;

    if (new_vert == mesh2.nPoint)
        mesh2.nPoint++;

    todelete[0] = e1;
    ntodelete = 1;
    addTria(v1, v2, new_vert, lab);

    chknadd(new_vert, v2);
    chknadd(v1, new_vert);
    for (nn = 0; nn < ntodelete; nn++)
        remFace(todelete[nn]);
    return 0;
} /*newTria*/

static void fill_eadj();
static void fill_tadj();
static int opt_func(int nfixed);

extern  char * ppMemory;
extern  PStrucFace2  *ptree2face;
extern  int StopAfterinitRegion;
extern int    nVRTglobal;
extern int    nTRIglobal;


/* error codes:
 *  0 - success
 * -1 - zero sized edge (error in user data)
 * -2 - internal search failed (most likely self intersection in front)
 ****************************************************************************/
int makeTria() {
    int i = 0, j, err = 0, smooth = 0;

    init();
    initRegion();

    if ( StopAfterinitRegion != 0 ) {
        free(ppMemory);
        free(ptree2face);
        return 0;
    }

    for (i = 1; i <= mesh2.nRegion; i++) {
        mesh2.nRPoint[i-1] = mesh2.nPoint;
        mesh2.nRTria[i-1] = mesh2.nTria;
        minrho = dist(tree2.face[0]->v1, tree2.face[0]->v2);

        for (j=0; j<tree2.nFace; j++) {
            if (minrho > dist(tree2.face[j]->v1, tree2.face[j]->v2))
                minrho = dist(tree2.face[j]->v1, tree2.face[j]->v2);
        }

        if (boolFAF)
            minrho *= 0.95;
        else
            minrho = 0.0;
        beta = 0.5;
        minrho *= beta;
        alpha = 2.0*minrho*(1.0-beta)/beta;
        alpha *= alpha;
        while (tree2.nFace > 0) {
            err = newTria(i);
            if (err)
                break;
            if (SHOWPROGRESS && (mesh2.nTria % 100 == 0)) {
                printf("Number of Point = %d Number of Tria = %d\n", mesh2.nPoint, mesh2.nTria);
                fflush(stdout);
            }
        }
        if (err)
            break;
        if (SHOWPROGRESS)
            printf("\n");
        fill_eadj();
        fill_tadj();

        smooth += opt_func(mesh2.nRPoint[i-1]);

        if (i != mesh2.nRegion) {
            initAddRegion(i + 1);
        }
        mesh2.nRPoint[i-1] = mesh2.nPoint - mesh2.nRPoint[i-1];
        mesh2.nRTria[i-1] = mesh2.nTria - mesh2.nRTria[i-1];
    }


    /*	printf("\nRESULT :  %5u     %5u    \n",mesh2.nPoint,mesh2.nTria);*/
    if (!smooth)
        regularity();

    test_quality();

    outMesh();
    if (!err) {
        free(ppMemory);
        free(ptree2face);
    } else {
        nVRTglobal = mesh2.nPoint;
        nTRIglobal = mesh2.nTria;
    }

    return err;
} /*makeTria*/

static void fill_eadj() {
    eadj.clear();
    eadj.resize(mesh2.nPoint);

    for (int i=0; i<mesh2.nTria; i++) {
        int a, b, c;
        a = mesh2.v1[i];
        b = mesh2.v2[i];
        c = mesh2.v3[i];

        eadj[a].insert(b);
        eadj[a].insert(c);
        eadj[b].insert(c);
        eadj[b].insert(a);
        eadj[c].insert(a);
        eadj[c].insert(b);
    }
}

static void fill_tadj() {
    tadj.clear();
    tadj.resize(mesh2.nPoint);

    for (int i=0; i<mesh2.nTria; i++) {
        int a;
        edge b;

        a    = mesh2.v1[i];
        b[0] = mesh2.v2[i];
        b[1] = mesh2.v3[i];

        tadj[a].insert(b);

        a    = mesh2.v2[i];
        b[0] = mesh2.v3[i];
        b[1] = mesh2.v1[i];
        tadj[a].insert(b);

        a    = mesh2.v3[i];
        b[0] = mesh2.v1[i];
        b[1] = mesh2.v2[i];
        tadj[a].insert(b);
    }
}

static int opt_func(int nfixed) {
    double x[2], z[2];

    std::vector<int> ps;

    double ds = 0.0;
    for (int p = nfixed; p < mesh2.nPoint; p++) {
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

    std::vector<vertex> dx(ps.size());
    std::vector<vertex> bkp(ps.size());

    for (size_t j = 0; j < ps.size(); j++) {
        int p = ps[j];
        bkp[j][0] = mesh2.x[p];
        bkp[j][1] = mesh2.y[p];
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
                func_xy(p, b, c, z, delta, 8);
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
            mesh2.x[p] += d*dx[j][0];
            mesh2.y[p] += d*dx[j][1];
        }
        if (SHOWPROGRESS) {
            printf(" Smoothing %d%%\n", s+1);
            fflush(stdout);
        }
    }
    if (SHOWPROGRESS)
        printf(" Smoothing done\n");

    int s = 0;
    for (int j = 0; j < mesh2.nTria; j++)
        if (func_q(mesh2.v1[j], mesh2.v2[j], mesh2.v3[j]) <= 0.0)
            s++;

    if (s) {
        fprintf(stderr,"aniAFT: Quality improvement failed, falling back to simple smoothing\n");
        for (size_t j = 0; j < ps.size(); j++) {
            int p = ps[j];
            mesh2.x[p] = bkp[j][0];
            mesh2.y[p] = bkp[j][1];
        }
    }

    return s;
}
