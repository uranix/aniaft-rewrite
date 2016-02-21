# CXX=clang++ -ferror-limit=5
CXX=g++ -fmax-errors=5
CXXFLAGS=-std=c++11 -O0 -g -Wall

CXXSRC=tree.cpp main.cpp mesh.cpp tria.cpp

CXXOBJ=${CXXSRC:.cpp=.o}

all: ${CXXOBJ} 
	c++ -o main $^ -lm

clean:
	rm -f *.o main
