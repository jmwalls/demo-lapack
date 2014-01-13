CC=g++
CFLAGS=-c -Wall
LDFLAGS=-llapack -lblas

all: lapack-demo

lapack-demo: demo.o
	$(CC) demo.o $(LDFLAGS) -o lapack-demo

demo.o: demo.cpp
	$(CC) $(CFLAGS) demo.cpp

clean:
	rm -rf *.o lapack-demo
