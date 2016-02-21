#ifndef __METRIC_H__
#define __METRIC_H__

#include <cmath>
#include <string>

struct coord {
    double x, y;
    coord() { }
    coord(double x, double y) : x(x), y(y) { }
    std::string str() const {
        return "{" + std::to_string(x) + ", " + std::to_string(y) + "}";
    }
};

struct Metric {
    static double distance(const coord &p1, const coord &p2) {
        return distance(p1.x, p1.y, p2.x, p2.y);
    }
    static double distance(double x1, double y1, double x2, double y2) {
        return std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    }
    virtual double size(const coord &p) const = 0;
};

struct UniformMetric : public Metric {
    const double h;
    UniformMetric(double h) : h(h) { }

    double size(const coord &) const { return h; }
};

#endif
