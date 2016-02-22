from aniaft import coord, Line, Bezier, DirectedSeg, Boundary, Triangulation

p = []
p.append(coord(0, 0))
p.append(coord(1, 0));
p.append(coord(1, 2));
p.append(coord(2, 1));
p.append(coord(4, 1.5));
p.append(coord(2, 0));
p.append(coord(3, -2));

print p
l1 = Line(1, 5, p)

s = []
s.append(Line(1, 5, p))
s.append(Line(3, 5, p))

reg0 = []
reg0.append(DirectedSeg(2, False))
reg0.append(DirectedSeg(1, False))
regs = [reg0]

bnd = Boundary(p, s, regs)
