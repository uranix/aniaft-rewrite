# CXX=clang++ -ferror-limit=500
CXX=g++ -fmax-errors=5
CXXFLAGS=-std=c++11 -O0 -g -Wall

SRC=aft2d.c memory.c refine.c region.c struct.c \
	tria.c user.c geometry.c main.c mesh.c \
	tree.c

OBJ=${SRC:.c=.o}

all: ${OBJ}
	c++ -o main ${OBJ} -lm -lgfortran

clean:
	rm -f *.o main
