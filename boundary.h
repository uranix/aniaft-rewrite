#ifndef H_BOUNDARY_MESH2D
#define H_BOUNDARY_MESH2D

#include "metric.h"

#include <cassert>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>

#ifndef SWIG

struct Transformation {
    const double m, sx, sy;
    Transformation(double m, double sx, double sy)
        : m(m), sx(sx), sy(sy)
    { }
    coord operator()(const coord &p) const {
        return coord(m * p.x + sx, m * p.y + sy);
    }
};

struct BoundingBox {
    bool initialized;
    double xmin, xmax, ymin, ymax;
    BoundingBox() : initialized(false) { }
    void add(const coord &o);
    Transformation toUnit() const;
    Transformation fromUnit() const;
};

#endif

struct Segment {
    int begCorner, endCorner;
    double tmin, tmax;
    int label;

    Segment(int begCorner, int endCorner, int label = 1, double tmin = 0, double tmax = 1)
        : begCorner(begCorner), endCorner(endCorner)
        , tmin(tmin), tmax(tmax), label(label)
    { }

    virtual coord operator()(const double t) const = 0;
    virtual std::string to_str() const {
        return "Segment("
            + std::to_string(begCorner) + ", "
            + std::to_string(endCorner) + ")";
    }
    virtual ~Segment() { }
};

struct Line : public Segment {
    const coord p0, p1;

    Line(int begCorner, int endCorner, const std::vector<coord> &corns, int label = 1)
        : Segment(begCorner, endCorner, label)
        , p0(corns[begCorner]), p1(corns[endCorner])
    { }

    coord operator()(const double t) const override {
        const double c0 = 1 - t;
        const double c1 = t;
        return coord(
                c0 * p0.x + c1 * p1.x,
                c0 * p0.y + c1 * p1.y
            );
    }

    std::string to_str() const override {
        return "Line[{"
            + p0.to_str() + ", "
            + p1.to_str() + "}]";
    }
};

struct Bezier : public Segment {
    const coord p0, p1, p2, p3;

    Bezier(
            int begCorner, int endCorner, const std::vector<coord> &corns,
            const coord &p1, const coord &p2, int label = 1
        )
        : Segment(begCorner, endCorner, label)
        , p0(corns[begCorner]), p1(p1), p2(p2), p3(corns[endCorner])
    { }

    coord operator()(const double t) const override {
        const double s = 1 - t;
        const double c0 =     s * s * s;
        const double c1 = 3 * t * s * s;
        const double c2 = 3 * t * t * s;
        const double c3 =     t * t * t;
        return coord(
                c0 * p0.x + c1 * p1.x + c2 * p2.x + c3 * p3.x,
                c0 * p0.y + c1 * p1.y + c2 * p2.y + c3 * p3.y
            );
    }

    std::string to_str() const override {
        return "BezierCurve[{"
            + p0.to_str() + ", "
            + p1.to_str() + ", "
            + p2.to_str() + ", "
            + p3.to_str() + "}]";
    }
};

struct DirectedSeg {
    int segment;
    bool reversed;

    DirectedSeg() { }

    DirectedSeg(int segment, bool reversed)
        : segment(segment), reversed(reversed)
    { }

    int vBeg(const std::vector<Segment *> &segs) const {
        if (!reversed)
            return segs[segment]->begCorner;
        else
            return segs[segment]->endCorner;
    }

    int vEnd(const std::vector<Segment *> &segs) const {
        if (reversed)
            return segs[segment]->begCorner;
        else
            return segs[segment]->endCorner;
    }
};

class Triangulation;

class Boundary {
    std::vector<coord> corners;
    std::vector<Segment *> segments;
    std::vector<std::vector<DirectedSeg> > regions;
    BoundingBox bb;
public:
    Boundary(
            const std::vector<coord> corners,
            const std::vector<Segment *> segments,
            const std::vector<std::vector<DirectedSeg> > regions
        )
        : corners(corners), segments(segments), regions(regions)
    {
        computeBoundingBox();
        checkRegions();
        checkCorners();
    }
private:
    void computeBoundingBox();
    bool checkCorner(int c, const coord &p) const;
    void checkCorners() const;
    void checkRegions() const;

    friend class Triangulation;
};

#endif
