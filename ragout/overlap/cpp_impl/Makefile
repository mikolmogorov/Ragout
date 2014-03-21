CC=g++
CFLAGS=-c -fPIC -Wall -std=c++0x -O3 -I${PYTHON_INCLUDE}
LDFLAGS=-shared
OUT_DIR=..

all: overlap_python

overlap_python: overlap.o fasta.o overlap_python.o
	$(CC) $(LDFLAGS) overlap.o fasta.o overlap_python.o -o ${OUT_DIR}/coverlap.so

overlap.o: overlap.cpp fasta.h overlap.h
	$(CC) $(CFLAGS) overlap.cpp -o overlap.o

overlap_python.o: overlap_python.cpp overlap.h
	$(CC) $(CFLAGS) overlap_python.cpp -o overlap_python.o 

fasta.o: fasta.h fasta.cpp
	$(CC) $(CFLAGS) fasta.cpp -o fasta.o

clean:
	rm *.o
	rm ${OUT_DIR}/coverlap.so
