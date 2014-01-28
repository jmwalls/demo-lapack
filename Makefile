CC=g++
CFLAGS=-c -Wall --std=gnu++0x
LDFLAGS=-llapack

INC=-I/usr/local/include/eigen3

all: lapack-demo

lapack-demo: demo.o linalg.o
	$(CC) demo.o linalg.o $(LDFLAGS) -o lapack-demo

demo.o: linalg.o demo.cpp
	$(CC) $(CFLAGS) $(INC) demo.cpp

linalg.o: linalg.h linalg.cpp
	$(CC) $(CFLAGS) $(INC) linalg.cpp

clean:
	rm -rf *.o lapack-demo
