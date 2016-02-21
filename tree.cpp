#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "tree.h"
#include "mesh.h"
#include "metric.h"

Tree::Tree() {
    root = new Node2d;
    root->entrycount = 0;
    root->parent = 0;
    root->nodelist[0] = 0;
    root->nodelist[1] = 0;
    root->nodelist[2] = 0;
    root->nodelist[3] = 0;
    root->firstentry = 0;
}

Tree::~Tree() {
    while (!faces.empty())
        remFace(faces.front());
    delete root;
}

int Tree::direction(double x, double y, double xc, double yc) {
    int dx=0, dy=0;
    if (x > xc)
        dx = 1;
    if (y > yc)
        dy = 2;
    return dx + dy;
}

void Tree::center(double *x, double *y, double side, int dir) {
    int dx=0, dy=0;

    if ((dir<0)||(dir>3))
        dir = 0;

    dx = dir%2;
    dy = dir/2;

    if (dx)
        x[0] += side;
    else
        x[0] -= side;

    if (dy)
        y[0] += side;
    else
        y[0] -= side;

    return;
}

void Tree::addFaceArray(PFace face) {
    int son,fath,i;
    PFace cface;

    faces.push_back(face);
    face->f = faces.size() - 1;
    son = faces.size();

aaa:
    fath=son/2;

    if ((fath>0) && (faces[son-1]->s < faces[fath-1]->s)) {
        cface = faces[son-1];
        faces[son-1]  = faces[fath-1];
        faces[fath-1] = cface;
        i = cface->f;
        cface->f = faces[son-1]->f;
        faces[son-1]->f = i;
        son=fath;
        if(son>1)
            goto aaa;
    }
}

void Tree::remFaceArray(int fath) {
    int son1,son2,i;
    PFace face;
    if((fath >= (int)faces.size()) || (fath < 0))  {
        fprintf(stderr, "aniAFT: delete face: tree is empty\n");
        return;
    }

    faces[fath] = faces.back();
    faces[fath]->f = fath;
    faces.pop_back();
aaa:i = 0;
    if (2*fath+1 < ((int)faces.size()-1)) {
        son1=2*fath+1;
        son2=son1++;
        if ((faces[fath]->s<=faces[son1]->s)&&(faces[fath]->s<=faces[son2]->s))
            i=0;
        if ((faces[fath]->s>=faces[son1]->s)&&(faces[son1]->s>=faces[son2]->s))
            i=son2;
        if ((faces[fath]->s>=faces[son2]->s)&&(faces[son2]->s>faces[son1]->s))
            i=son1;
        if ((faces[son1]->s>=faces[fath]->s)&&(faces[fath]->s>faces[son2]->s))
            i=son2;
        if ((faces[son2]->s>=faces[fath]->s)&&(faces[fath]->s>faces[son1]->s))
            i=son1;
    }
    if(2*fath+1 == ((int)faces.size()-1)){
        son1=2*fath+1;
        if ((faces[fath]->s>faces[son1]->s))
            i=son1;
        else
            i=0;
    }

    if(i)
    {
        face=faces[i];
        faces[i]=faces[fath];
        faces[fath]=face;
        son1=face->f;
        face->f=faces[i]->f;
        faces[i]->f=son1;
        fath=i;
        goto aaa;
    }
}

PFace Tree::addFace(const Mesh &mesh, int v1, int v2) {
    int         i,d;
    double      x1=-1.0,y1=-1.0;
    double      xc=0.5,yc=0.5,side=0.5;
    PFace face;
    PNode2d old_node,node;
    entrychain  *faceentry, *new_faceentry;

    node = root;
    face = new Face;
    face->v1 = v1;
    face->v2 = v2;

    const Vertex &p1 = mesh.pts[v1];
    const Vertex &p2 = mesh.pts[v2];

    double x = (p1.x + p2.x)/2.;
    double y = (p1.y + p2.y)/2.;

    face->x=x;
    face->y=y;
    face->s = Metric::distance(p1.x, p1.y, p2.x, p2.y);

    if ((x<0.)||(x>1.)||(y<0.)||(y>1)) {
        printf("x= %lf,  y= %lf ", x, y);
        fprintf(stderr, "aniAFT: out of bouding box (wrong orientation in input?)\n");
    }
    while (side>face->s) {
        d = direction(x,y,xc,yc);
        side *= 0.5;
        center(&xc,&yc,side,d);
        old_node = node;
        node = node->nodelist[d];
        if (!node) {
            node = new Node2d;
            for (i=0;i<4;i++)
                node->nodelist[i] = 0;
            node->parent = old_node;
            node->entrycount = 0;
            node->firstentry = 0;
            old_node->nodelist[d] = node;
        }
    }
    faceentry = node->firstentry;
    for (i=0;i<node->entrycount;i++) {
        x1 = faceentry->face->x;
        y1 = faceentry->face->y;
        if ((x1==x)&&(y1==y)&&(v1==faceentry->face->v1)&&(v2==faceentry->face->v2)) {
            printf("x,  y:  %lf, %lf\nx1, y1: %lf, %lf\nv1:   %d, v2:   %d\n v1_f: %d, v2_f: %d\n", x, y, x1, y1, v1, v2, faceentry->face->v1, faceentry->face->v2);
            fprintf(stderr, "aniAFT: dups in front (corrupted front?)\n");
        }
        faceentry = faceentry->next;
    }
    faceentry = node->firstentry;
    new_faceentry = new entrychain;
    node->firstentry = new_faceentry;
    new_faceentry->face = face;
    new_faceentry->next = faceentry;
    node->entrycount++;
    addFaceArray(face);

    return face;
}

void Tree::remFace(PFace face) {
    int         i,d;
    double      x,y;
    double      xc=0.5,yc=0.5,side=0.5;
    PNode2d node, old_node;
    entrychain  *faceentry, *old_faceentry;

    remFaceArray(face->f);
    x = face->x;
    y = face->y;
    node = root;
    while (side>face->s) {
        d = direction(x,y,xc,yc);
        side *= 0.5;
        center(&xc,&yc,side,d);
        node = node->nodelist[d];
        if (!node) {
            fprintf(stderr, "aniAFT: delete face: no such element\n");
            return;
        }
    }
    faceentry = node->firstentry;
    if (faceentry->face == face) {
        node->firstentry = faceentry->next;
        delete face;
        delete faceentry;
        node->entrycount--;
    } else {
        do {
            old_faceentry = faceentry;
            faceentry = faceentry->next;
        } while ((faceentry)&&(faceentry->face != face));
        if (faceentry) {
            old_faceentry->next = faceentry->next;
            delete face;
            delete faceentry;
            node->entrycount--;
        } else {
            fprintf(stderr, "aniAFT: delete face: no such element\n");
            return;
        }
    }
    while ((node->entrycount==0)&&(!node->nodelist[0]&&!node->nodelist[1]&&
                !node->nodelist[2]&&!node->nodelist[3])&&(node!=root)) {
        old_node = node;
        node = node->parent;
        for (i=0;i<4;i++)
            if (node->nodelist[i]==old_node) break;
        node->nodelist[i] = 0;
        delete old_node;
    }
    return;
} /* remFace */

int Tree::sectQuad() const {
    /*  need  Quad     : tree.xc;  tree.yc;
        tree.side;
        Vicinity : tree.xVicinity;  tree.yVicinity;
        tree.sVicinity;
        */
    if(xc+3*side < xVicinity-sVicinity) return 0;
    if(xc-3*side > xVicinity+sVicinity) return 0;
    if(yc+3*side < yVicinity-sVicinity) return 0;
    if(yc-3*side > yVicinity+sVicinity) return 0;

    return 1;
}

void Tree::vicinityFacesRec(PNode2d node) {
    /*  NEED    xc=0.5; yc=0.5; side=0.5;
        vicinityFaces empty
        sVicinity  xVicinity yVicinity
        */
    int         i;
    double      x,y,s;
    entrychain  *faceentry;

    if (!node)
        return;

    for (faceentry = node->firstentry;faceentry;faceentry=faceentry->next) {
        x = faceentry->face->x;
        y = faceentry->face->y;
        s = Metric::distance(x,y,xVicinity,yVicinity)-0.5*faceentry->face->s;
        if (s<=sVicinity)
            vicinityFaces.push_back(faceentry->face);
    }
    side *= 0.5;
    x = xc;
    y = yc;
    s = side;
    if (node->nodelist[0]) {
        center(&xc, &yc, side, 0);
        if (sectQuad())
            vicinityFacesRec(node->nodelist[0]);
    }
    for (i=1;i<4;i++) {
        if (!node->nodelist[i])
            continue;
        side = s;
        xc = x;
        yc = y;
        center(&xc, &yc, side, i);
        if (sectQuad())
            vicinityFacesRec(node->nodelist[i]);
    }
}

void Tree::buildVicinityFaces(double x, double y, double size) {
    xc=0.5;
    yc=0.5;
    side=0.5;

    vicinityFaces.clear();

    sVicinity = size;
    xVicinity = x;
    yVicinity = y;

    vicinityFacesRec(root);
}
