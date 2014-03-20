CC=g++
CFLAGS=-c -Wall -std=c++0x -O3
LDFLAGS=

BIN_DIR=../bin

all: overlap

overlap: overlap.o fasta.o
	$(CC) $(LDFLAGS) overlap.o fasta.o -o ${BIN_DIR}/overlap

overlap.o: overlap.cpp fasta.h
	$(CC) $(CFLAGS) -c overlap.cpp -o overlap.o

fasta.o: fasta.h fasta.cpp
	$(CC) $(CFLAGS) -c fasta.cpp -o fasta.o

clean:
	rm *.o
	rm ${BIN_DIR}/overlap
