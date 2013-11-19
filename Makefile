#By default the compiler is G++. Though GCC may be more efficient.
CC=g++
CFLAGS=-c -Wall -frounding-math -g -O3 #-O3 
LDFLAGS= -lCGAL -llapack -lblas -lmpfr -lgmp
SRC = $(wildcard *.cc)
EXECUTABLE = pcdqb

all: $(SRC:%.cc=%.o)
	$(CC) $(SRC:%.cc=%.o) $(LDFLAGS) -o $(EXECUTABLE)

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)


