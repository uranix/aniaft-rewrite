#ifndef H_BOUNDARY_MESH2D
#define H_BOUNDARY_MESH2D

#include "metric.h"

#include <cassert>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>

struct Transformation {
    const double m, sx, sy;
    Transformation(double m, double sx, double sy)
        : m(m), sx(sx), sy(sy)
    { }
    coord operator()(const coord &p) const {
        return coord(m * p.x + sx, m * p.y + sy);
    }
};

struct Segment {
    int begCorner, endCorner;
    double tmin, tmax;
    int label;

    Segment(int begCorner, int endCorner, int label = 1, double tmin = 0, double tmax = 1)
        : begCorner(begCorner), endCorner(endCorner)
        , tmin(tmin), tmax(tmax)
    { }

    virtual coord operator()(const double t) const = 0;
    virtual std::string str() const {
        return "Segment("
            + std::to_string(begCorner) + ", "
            + std::to_string(endCorner) + ")";
    }
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

    std::string str() const override {
        return "Line[{"
            + p0.str() + ", "
            + p1.str() + "}]";
//            + std::to_string(begCorner) + ":" + std::to_string(endCorner) + ")";
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

    std::string str() const override {
        return "BezierCurve[{"
            + p0.str() + ", "
            + p1.str() + ", "
            + p2.str() + ", "
            + p3.str() + "}]";
//            + std::to_string(begCorner) + ":" + std::to_string(endCorner) + ")";
    }
};

struct BoundingBox {
    bool initialized;
    double xmin, xmax, ymin, ymax;
    BoundingBox() : initialized(false) { }
    void add(const coord &p) {
        if (!initialized) {
            initialized = true;
            xmin = xmax = p.x;
            ymin = ymax = p.y;
            return;
        }
        if (p.x < xmin)
            xmin = p.x;
        if (p.x > xmax)
            xmax = p.x;
        if (p.y < ymin)
            ymin = p.y;
        if (p.y > ymax)
            ymax = p.y;
    }
    Transformation toUnit() const {
        double xlft = 1.1 * xmin - 0.1 * xmax;
        double ybot = 1.1 * ymin - 0.1 * ymax;
        double size = 1.2 * std::max(xmax - xmin, ymax - ymin);
        return Transformation(
                1. / size, -xlft / size, -ybot / size
            );
    }
    Transformation fromUnit() const {
        double xlft = 1.1 * xmin - 0.1 * xmax;
        double ybot = 1.1 * ymin - 0.1 * ymax;
        double size = 1.2 * std::max(xmax - xmin, ymax - ymin);
        return Transformation(size, xlft, ybot);
    }
};

struct Link {
    int segment;
    bool reversed;

    Link(int segment, bool reversed)
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

struct Boundary {
    std::vector<coord> corners;
    std::vector<Segment *> segments;
    std::vector<std::vector<Link> > regions;

    BoundingBox bb;

    Boundary(
            const std::vector<coord> corners,
            const std::vector<Segment *> segments,
            const std::vector<std::vector<Link> > regions
        )
        : corners(corners), segments(segments), regions(regions)
    {
        computeBoundingBox();
        checkRegions();
        checkCorners();
    }

    void computeBoundingBox() {
        for (const Segment *s : segments) {
            double dt = (s->tmax - s->tmin) / 100;
            for (double t = s->tmin; t < s->tmax; t += dt)
                bb.add((*s)(t));
        }
    }

    bool checkCorner(int c, const coord &p) const {
        double tol = 1e-8;
        const auto transf = bb.toUnit();
        const coord &p1 = transf(corners[c]);
        const coord &p2 = transf(p);

        return Metric::distance(p1.x, p1.y, p2.x, p2.y) < tol;
    }

    void checkCorners() const {
        assert(bb.initialized);
        for (const Segment *s : segments) {
            bool ret;
            ret = checkCorner(s->begCorner, (*s)(s->tmin));
            if (!ret)
                throw std::logic_error("The beginning of " + s->str() + " mismatch the corner point " + corners[s->begCorner].str());
            ret = checkCorner(s->endCorner, (*s)(s->tmax));
            if (!ret)
                throw std::logic_error("The end of " + s->str() + " mismatch the corner point " + corners[s->endCorner].str());
        }
    }

    void checkRegions() const {
        for (size_t k = 0; k < regions.size(); k++) {
            const auto &r = regions[k];
            bool newLoop = true;
            int s, n;
            for (size_t j = 0; j < r.size(); j++) {
                if (newLoop) {
                    s = r[j].vBeg(segments);
                    n = r[j].vEnd(segments);
                    newLoop = false;
                } else {
                    if (r[j].vBeg(segments) != n)
                        throw std::logic_error("Inspect segments " + std::to_string(j-1) + " and " + std::to_string(j)
                                + " in region " + std::to_string(k));
                    n = r[j].vEnd(segments);
                    if (n == s)
                        newLoop = true;
                }
            }
            if (n != s)
                throw std::logic_error("Region " + std::to_string(k) + " is not closed");
        }
    }
};

#endif
