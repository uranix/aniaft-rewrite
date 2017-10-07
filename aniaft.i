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

%feature("kwargs") Segment;

%include "metric.h"

%include "mesh.h"

%ignore Mesh;

%template(Coords) std::vector<coord>;
%template(Segments) std::vector<Segment *>;
%template(Chain) std::vector<DirectedSeg>;
%template(Regions) std::vector<std::vector<DirectedSeg> >;
%template(Triangles) std::vector<Triangle>;
%template(Edges) std::vector<Edge>;

%feature("director") Segment;

%include "boundary.h"

%include "tria.h"
