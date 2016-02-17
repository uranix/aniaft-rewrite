# CXX=clang++ -ferror-limit=500
CXX=g++ -fmax-errors=5
CXXFLAGS=-std=c++11 -O3 -g -Wall

F77=gfortran
F77FLAGS=-O3 -g -Wall

CXXSRC=aft2d.cpp memory.cpp refine.cpp region.cpp mesh.cpp \
	tria.cpp user.cpp geometry.cpp mesh.cpp \
	tree.cpp

F77SRC=main.f graph.f

CXXOBJ=${CXXSRC:.cpp=.o}
F77OBJ=${F77SRC:.f=.o}

all: ${F77OBJ} ${CXXOBJ} 
	c++ -o main $^ -lm -lgfortran

clean:
	rm -f *.o main
