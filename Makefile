CXX=g++
CXXFLAGS=-std=c++11 -O3 -Wall

SRC=aft2d.c memory2.c refine2.c region2.c struct2.c \
	tria2.c user2.c user.c main.c mesh.c \
	tree2.c

OBJ=${SRC:.c=.o}

all: ${OBJ}
	c++ -o main ${OBJ} -lm -lgfortran
