#include <stdlib.h>
#include "tree2.h"
#include "region2.h"
#include "refine2.h"

#include "tree2.h"
#include "region2.h"
#include "refine2.h"

#include <cassert>
#include <map>
#include <algorithm>

extern  StrucMesh2  mesh;
extern  int  boolRegul;
extern  int  boolFAF;

double angle( int v1, int v2, int v ) {
    double x1,x2,y1,y2,xx,yy;

    x1=mesh.pts[v1].x;
    y1=mesh.pts[v1].y;
    x2=mesh.pts[v2].x;
    y2=mesh.pts[v2].y;

    xx=mesh.pts[v].x;
    yy=mesh.pts[v].y;

    x2-=xx;
    x1-=xx;
    y1-=yy;
    y2-=yy;

    double p = x1*x2+y1*y2;
    p /= sqrt((x1*x1+y1*y1) * (x2*x2+y2*y2));

    return acos(p);
}

void test_quality()
{
    double da=M_PI/6,min=10.0,max=0.0;
    double min_angle=3.0,max_angle=0.0;

    std::map<size_t, int> n;

    for (const auto &tr : mesh.tri) {
        int p1=tr.v1;
        int p2=tr.v2;
        int p3=tr.v3;

        double a1=angle(p2,p3,p1);
        double a2=angle(p3,p1,p2);
        double a3=angle(p1,p2,p3);

        for(int j = 0; j < 6; j++){
            if ((j*da<a1) && (a1<=j*da+da))
                n[j]++;
            if ((j*da<a2) && (a2<=j*da+da))
                n[j]++;
            if ((j*da<a3) && (a3<=j*da+da))
                n[j]++;
            if (a1>max_angle)
                max_angle=a1;
            if (a1<min_angle)
                min_angle=a1;
            if (a2>max_angle)
                max_angle=a2;
            if (a2<min_angle)
                min_angle=a2;
            if (a3>max_angle)
                max_angle=a3;
            if (a3<min_angle)
                min_angle=a3;
        }
        double x1=mesh.pts[p1].x;
        double y1=mesh.pts[p1].y;
        double x2=mesh.pts[p2].x;
        double y2=mesh.pts[p2].y;
        double x3=mesh.pts[p3].x;
        double y3=mesh.pts[p3].y;

        double xc=0.3333333*(x1+x2+x3);
        double yc=0.3333333*(y1+y2+y3);

        double size, pp;

        if (!boolFAF)
            size=1.0/sizeFace(xc,yc);
        else
            size=3.0/(distance(x1,y1,x2,y2)+distance(x2,y2,x3,y3)+distance(x3,y3,x1,y1));

        pp=size*distance(x1,y1,x2,y2);
        if(pp<min)
            min=pp;
        if(pp>max)
            max=pp;
        pp=size*distance(x3,y3,x2,y2);
        if(pp<min)
            min=pp;
        if(pp>max)
            max=pp;
        pp=size*distance(x1,y1,x3,y3);
        if(pp<min)
            min=pp;
        if(pp>max)
            max=pp;
    }

    for (const auto p : mesh.pts)
        n[p.neib.size()]++;

    printf("Neigbor number for nP = %7lu nT = %7lu\n", mesh.pts.size(), mesh.tri.size());
    for (const auto &kv : n) {
        int i = kv.first;
        int ni = kv.second;
        printf("neig =%6d  nPoint =%6d  perc = %lf  perc = %lf\n", i, ni,
                (double)ni/mesh.pts.size(),(double)ni/(mesh.pts.size()-n[0]));
    }
}

void smoothing() {
    for (Point &p : mesh.pts) {
        double xx=0.;
        double yy=0.;
        for (const auto &nei : p.neib) {
            double x = mesh.pts[nei].x;
            double y = mesh.pts[nei].y;
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

void delPoint(int n) {
    mesh.pts[n].remove = true;
}


void changePoint(int v1, int v2) {
    mesh.pts[v1].x = 0.5*(mesh.pts[v1].x + mesh.pts[v2].x);
    mesh.pts[v1].y = 0.5*(mesh.pts[v1].y + mesh.pts[v2].y);
}

void delTria(int n) {
    mesh.tri[n].remove = true;
}

void changeTria(int iTria, int v1, int v2, int v3) {
    mesh.tri[iTria].v1 = v1;
    mesh.tri[iTria].v2 = v2;
    mesh.tri[iTria].v3 = v3;
}

void sortNeigbor(int iVert) {
    std::vector<std::pair<int, double> > na;

    double x0 = mesh.pts[iVert].x;
    double y0 = mesh.pts[iVert].y;

    for (const int nei : mesh.pts[iVert].neib) {
        double x = mesh.pts[nei].x - x0;
        double y = mesh.pts[nei].y - y0;
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

    std::transform(na.begin(), na.end(), mesh.pts[iVert].neib.begin(),
            [] (const std::pair<int, double> &v) { return v.first; });
}

void calcNeigTria() {
    for (Point &p : mesh.pts)
        p.neibTria.clear();
    int i = 0;
    for (const auto tr : mesh.tri) {
        mesh.pts[tr.v1].neibTria.push_back(i);
        mesh.pts[tr.v2].neibTria.push_back(i);
        mesh.pts[tr.v3].neibTria.push_back(i);
        i++;
    }
}

void calcNeigbor() {
    std::vector<int> vert;

    for (Point &p : mesh.pts) {
        if (!p.skip_neib)
            p.neib.clear();
    }

    for (int i = 0; i < (int)mesh.pts.size(); i++) {
        Point &p = mesh.pts[i];
        if (p.skip_neib)
            continue;
        vert.clear();
        for (const int itr : p.neibTria) {
            const Triangle &tr = mesh.tri[itr];
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

void calcEdge() {
    for(int i = 0; i < (int)mesh.pts.size(); i++) {
        Point &p = mesh.pts[i];
        if (p.skip_neib)
            continue;
        for (const int vert : p.neib) {
            if (vert <= i)
                continue;
            if (mesh.pts[vert].skip_neib)
                continue;

            mesh.vb.push_back(i);
            mesh.ve.push_back(vert);
            bool swap = true;
            for (size_t itr : p.neibTria) {
                const Triangle &tr = mesh.tri[itr];
                if (tr.v1 == vert || tr.v2 == vert || tr.v3 == vert) {
                    if (swap) {
                        mesh.tria1.push_back(itr);
                        swap = false;
                    } else
                        mesh.tria2.push_back(itr);
                }
            }
            swap = true;
            int itr = mesh.tria1.back();
            const Triangle &tr = mesh.tri[itr];
            if (tr.v1 == i && tr.v2 == vert)
                swap = false;
            else if (tr.v2 == i && tr.v3 == vert)
                swap = false;
            else if (tr.v3 == i && tr.v1 == vert)
                swap = false;
            if (swap)
                std::swap(mesh.tria1.back(), mesh.tria2.back());
        }
    }
}

void initRegularity() {
    calcNeigTria();
}

template<class Cont, class Pred>
void drop(Cont &cont, Pred pred) {
    cont.erase(std::remove_if(cont.begin(), cont.end(), pred), cont.end());
}

void pack() {
    drop(mesh.tri, [] (const Triangle &t) { return t.remove; });

    std::vector<int> newidx;

    int idx = 0;
    for (const Point &p : mesh.pts) {
        if (p.remove)
            newidx.push_back(-1);
        else {
            newidx.push_back(idx);
            idx++;
        }
    }

    for (auto &tr : mesh.tri) {
        tr.v1 = newidx[tr.v1];
        tr.v2 = newidx[tr.v2];
        tr.v3 = newidx[tr.v3];
        assert(tr.v1 >= 0 && tr.v2 >= 0 && tr.v3 >= 0);
    }

    drop(mesh.pts, [] (const Point &p) { return p.remove; });
}

// Deletes points with 3 or 4 neighbors
void deletePoint() {
    for(int i = 0; i < (int)mesh.pts.size(); i++) {
        Point &p = mesh.pts[i];
        int j = p.neib.size();
        if (j == 3) {
            delPoint(i);
            sortNeigbor(i);
            changeTria(p.neibTria[0], p.neib[1], p.neib[0], p.neib[2]);
            delTria   (p.neibTria[1]);
            delTria   (p.neibTria[2]);
        }
        else if( j == 4 ){
            delPoint(i);
            sortNeigbor(i);
            int n1 = mesh.pts[p.neib[0]].neib.size() - 7;
            int n2 = mesh.pts[p.neib[1]].neib.size() - 7;
            int n3 = mesh.pts[p.neib[2]].neib.size() - 7;
            int n4 = mesh.pts[p.neib[3]].neib.size() - 7;
            int sum  = (n1+1)*(n1+1) + (n3+1)*(n3+1) + n2*n2 + n4*n4;
            int swap = (n2+1)*(n2+1) + (n4+1)*(n4+1) + n1*n1 + n3*n3;
            /*printf("swap %4d%4d n(%4d%4d%4d%4d) \n",sum,swap,n1,n3,n2,n4);*/
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

    return;
}

// Merges two points with 5 neibs
void mergePoint() {
    calcEdge();
    for(size_t i = 0; i < mesh.vb.size(); i++){
        int v1 = mesh.vb[i];
        int v2 = mesh.ve[i];
        int n1 = mesh.pts[v1].neib.size();
        int n2 = mesh.pts[v2].neib.size();
        if( n1 == 5 && n2 == 5 ){
            // Dont merge futher these points
            mesh.pts[v1].neib.push_back(-1);
            mesh.pts[v2].neib.push_back(-1);
            changePoint(v1, v2);
            delPoint(v2);
            for (const int itr : mesh.pts[v2].neibTria) {
                Triangle &tr = mesh.tri[itr];
                if( tr.v1 == v2 )
                    tr.v1 = v1;
                if( tr.v2 == v2 )
                    tr.v2 = v1;
                if( tr.v3 == v2 )
                    tr.v3 = v1;
            }
            delTria(mesh.tria1[i]);
            delTria(mesh.tria2[i]);
        }
    }

    pack();
    calcNeigTria();
    calcNeigbor();
}

void splitPoint() {
    for (int i=0; i < (int)mesh.pts.size(); i++){
        Point &p = mesh.pts[i];
        int n1 = p.neib.size();
        if( n1 <= 7 )
            continue;
        int n2 = n1/2;
        n1 = n1 - n2;
        sortNeigbor(i);

        double size;

        if (!boolFAF)
            size = sizeFace(p.x, p.y);
        else {
            const Point &p2 = mesh.pts[p.neib[n2]];
            size = distance(p.x, p.y, p2.x, p2.y);
        }

        vert = mesh.nPoint;
        mesh.neigbor[0][vert] = 1;
        addPoint(mesh.x[i]+0.1*size,mesh.y[i]);

        for(j=1;j<=n1+n2;j++){
            if( mesh.neigbor[0][ mesh.neigbor[j][i] ] != -1 )
                mesh.neigbor[0][ mesh.neigbor[j][i] ] = 1;
        }
        lab=-1;
        for(j=1;j<n1;j++){
            it=mesh.neigTria[j][i]; if (lab==-1) lab=mesh.label[it];
            if (lab != mesh.label[it])  fprintf(stderr, "aniAFT: different materials inside region\n");
            changeTria(mesh.neigTria[j][i],
                    mesh.neigbor[j+1][i],mesh.neigbor[j][i],i);
        }
        for(j=n1+1;j<n1+n2;j++){
            it=mesh.neigTria[j][i]; if (lab==-1) lab=mesh.label[it]; 
            if (lab != mesh.label[it])  fprintf(stderr, "aniAFT: different materials inside region\n");
            changeTria(mesh.neigTria[j][i],
                    mesh.neigbor[j+1][i],mesh.neigbor[j][i],vert);
        }
        it=mesh.neigTria[n1][i]; if (lab==-1) lab=mesh.label[it]; 
        if (lab != mesh.label[it])  fprintf(stderr, "aniAFT: different materials inside region\n");
        changeTria(mesh.neigTria[n1][i],
                mesh.neigbor[n1][i],i,vert);
        it=mesh.neigTria[n1+n2][i]; if (lab==-1) lab=mesh.label[it]; 
        if (lab != mesh.label[it])  fprintf(stderr, "aniAFT: different materials inside region\n");
        changeTria(mesh.neigTria[n1+n2][i],
                mesh.neigbor[n1+1][i],mesh.neigbor[n1][i],vert);
        addTria(mesh.neigbor[n1+n2][i],vert,i,lab);
        addTria(mesh.neigbor[n1+n2][i],i,mesh.neigbor[1][i],lab);
    }
    /*   pack();*/
    calcNeigTria();
    calcNeigbor();
    for(i=0;i<5;i++)
        smoothing();

    return;
}/*splitPoint*/


void swapEdge()
{
    int  i,vb,ve,v1,v2,nb,ne,n1,n2,tria1,tria2,sum,swap;

    calcEdge();
    for(i=0;i<mesh.nEdge;i++){
        vb = mesh.vb[i];
        ve = mesh.ve[i];
        tria1 = mesh.tria1[i];
        tria2 = mesh.tria2[i];
        v1 = mesh.v1[tria1];
        if( v1 == vb || v1 == ve ){
            v1 = mesh.v2[tria1];
            if( v1 == vb || v1 == ve ){
                v1 = mesh.v3[tria1];
            }
        }
        v2 = mesh.v1[tria2];
        if( v2 == vb || v2 == ve ){
            v2 = mesh.v2[tria2];
            if( v2 == vb || v2 == ve ){
                v2 = mesh.v3[tria2];
            }
        }
        if( mesh.neigbor[0][v1] == -1 || mesh.neigbor[0][v2] == -1 )
            continue;
        if( mesh.neigTria[0][vb] == -1 && mesh.neigTria[0][ve] == -1 )
            continue;
        nb = mesh.neigbor[0][vb] - 6;
        ne = mesh.neigbor[0][ve] - 6;
        n1 = mesh.neigbor[0][v1] - 6;
        n2 = mesh.neigbor[0][v2] - 6;
        if( nb < 1 && ne < 1 && n1 < 1 && n2 < 1 )
            continue;
        sum = nb*nb + ne*ne + n1*n1 + n2*n2;
        swap = (nb-1)*(nb-1) + (ne-1)*(ne-1) + (n1+1)*(n1+1) + (n2+1)*(n2+1);
        /*printf("swap %4d%4d v(%4d%4d%4d%4d) n(%4d%4d%4d%4d) \n",sum,swap,vb,ve,v1,v2,nb,ne,n1,n2);*/
        if( swap < sum ){
            mesh.neigbor[0][vb] = nb + 5;
            mesh.neigbor[0][ve] = ne + 5;
            mesh.neigbor[0][v1] = n1 + 7;
            mesh.neigbor[0][v2] = n2 + 7;
            mesh.neigTria[0][vb] = -1;
            mesh.neigTria[0][ve] = -1;
            mesh.neigTria[0][v1] = -1;
            mesh.neigTria[0][v2] = -1;
            changeTria(tria1,v2,v1,vb);
            changeTria(tria2,v1,v2,ve);
        }
    }
    free(mesh.vb);
    free(mesh.ve);
    free(mesh.tria1);
    free(mesh.tria2);

    calcNeigTria();
    calcNeigbor();

    return;
}/*swapEdge*/


void regularity()
{
    int  i;

    if( boolRegul ){

        initRegularity();
        for(i=0;i<5;i++)
            smoothing();
        /*test_quality();
          getch();*/

        deletePoint();
        for(i=0;i<5;i++)
            smoothing();
        /*test_quality();
          getch();*/

        mergePoint();
        mergePoint();
        mergePoint();
        for(i=0;i<5;i++)
            smoothing();
        /*test_quality();
          getch();*/

        splitPoint();
        splitPoint();
        splitPoint();
        /*test_quality();
          printf("SWAP:\n");
          getch();*/

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

        /*test_quality();
          getch();*/

    }/*if( boolRegul )*/

    for(i=0;i<mesh.nSmooth;i++)
        smoothing();

    return;
}/*regularity*/
