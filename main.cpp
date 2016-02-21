#include "tria.h"

#include <iostream>

struct MyMetric : public Metric {
    double size(const coord &c) const override {
        return 0.03 * (3 + sin(4 * c.x) + sin(5 * c.y) );
    }
};

int main() {
    std::vector<coord> p;

    p.emplace_back(0, 0);
    p.emplace_back(1, 0);
    p.emplace_back(1, 2);
    p.emplace_back(2, 1);
    p.emplace_back(4, 1.5);
    p.emplace_back(2, 0);
    p.emplace_back(3, -2);

    std::vector<Segment *> s;

    s.push_back(new Line(1, 5, p)); // 0
    s.push_back(new Line(3, 5, p)); // 1
    s.push_back(new Line(2, 3, p)); // 2
    s.push_back(new Line(5, 6, p)); // 3
    s.push_back(new Line(3, 4, p)); // 4
    s.push_back(new Bezier(              // 5
                2, 6, p,
                coord(-2, 1), coord(-2, -1)
            ));
    s.push_back(new Line(6, 4, p)); // 6
    s.push_back(new Line(4, 2, p)); // 7
    s.push_back(new Bezier(              // 8
                1, 0, p,
                coord(1, .667), coord(0, .667)
            ));
    s.push_back(new Bezier(              // 9
                0, 1, p,
                coord(0, -.667), coord(1, -.667)
            ));

    std::vector<std::vector<Link> > regs(3);

    regs[0].emplace_back(2, false);
    regs[0].emplace_back(1, false);
    regs[0].emplace_back(3, false);
    regs[0].emplace_back(5, true);
    regs[0].emplace_back(8, false);
    regs[0].emplace_back(9, false);

    regs[1].emplace_back(2, true);
    regs[1].emplace_back(7, true);
    regs[1].emplace_back(4, true);

    regs[2].emplace_back(4, false);
    regs[2].emplace_back(6, true);
    regs[2].emplace_back(3, true);
    regs[2].emplace_back(1, true);

    Boundary bnd(p, s, regs);
    Triangulation tr(bnd, MyMetric());

    tr.generate();
    tr.mesh.saveVtk("out.vtk");

    return 0;
}
