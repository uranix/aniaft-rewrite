# CXX=clang++ -ferror-limit=500
CXX=g++ -fmax-errors=5
CXXFLAGS=-std=c++11 -O3 -g -Wall

F77=gcc
F77FLAGS=-O3 -g -Wall

CXXSRC=tree.cpp aft2d.cpp refine.cpp mesh.cpp \
	tria.cpp user.cpp mesh.cpp region.cpp

CSRC=main_boundary_wing.c demo.c crv_model.c

CXXOBJ=${CXXSRC:.cpp=.o}
COBJ=${CSRC:.c=.o}

all: ${COBJ} ${CXXOBJ} 
	c++ -o main $^ -lm -lf2c

clean:
	rm -f *.o main
