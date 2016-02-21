#include "boundary.h"

void BoundingBox::add(const coord &p) {
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

Transformation BoundingBox::toUnit() const {
    double xlft = 1.1 * xmin - 0.1 * xmax;
    double ybot = 1.1 * ymin - 0.1 * ymax;
    double size = 1.2 * std::max(xmax - xmin, ymax - ymin);
    return Transformation(
            1. / size, -xlft / size, -ybot / size
        );
}
Transformation BoundingBox::fromUnit() const {
    double xlft = 1.1 * xmin - 0.1 * xmax;
    double ybot = 1.1 * ymin - 0.1 * ymax;
    double size = 1.2 * std::max(xmax - xmin, ymax - ymin);
    return Transformation(size, xlft, ybot);
}

void Boundary::computeBoundingBox() {
    for (const Segment *s : segments) {
        double dt = (s->tmax - s->tmin) / 100;
        for (double t = s->tmin; t < s->tmax; t += dt)
            bb.add((*s)(t));
    }
}

bool Boundary::checkCorner(int c, const coord &p) const {
    double tol = 1e-8;
    const auto transf = bb.toUnit();
    const coord &p1 = transf(corners[c]);
    const coord &p2 = transf(p);

    return Metric::distance(p1.x, p1.y, p2.x, p2.y) < tol;
}

void Boundary::checkCorners() const {
    assert(bb.initialized);
    for (const Segment *s : segments) {
        bool ret;
        ret = checkCorner(s->begCorner, (*s)(s->tmin));
        if (!ret)
            throw std::logic_error("The beginning of " + s->to_str() + " mismatch the corner point " + corners[s->begCorner].to_str());
        ret = checkCorner(s->endCorner, (*s)(s->tmax));
        if (!ret)
            throw std::logic_error("The end of " + s->to_str() + " mismatch the corner point " + corners[s->endCorner].to_str());
    }
}

void Boundary::checkRegions() const {
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
