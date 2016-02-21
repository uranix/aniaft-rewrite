/* vim: set ft=cpp: */
%module(directors="1") aniaft

%{
#include "metric.h"
#include "boundary.h"
#include "tria.h"
#include "mesh.h"
%}

%include "std_string.i"
%include "std_vector.i"

%feature("director") Metric;

%include "metric.h"

%include "mesh.h"

%ignore Mesh;

%template(Coords) std::vector<coord>;
%template(Segments) std::vector<Segment *>;
%template(Chain) std::vector<Link>;
%template(Regions) std::vector<std::vector<Link> >;
%template(Triangles) std::vector<Triangle>;
%template(Edges) std::vector<Edge>;

%feature("director") Segment;

%include "boundary.h"

%ignore edge::operator[];

%include "tria.h"
