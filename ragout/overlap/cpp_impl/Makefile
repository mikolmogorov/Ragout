CFLAGS=-c -fPIC -Wall -O3 -I${PYTHON_INCLUDE}
LDFLAGS=-shared
OUT_DIR=..

.PHONY: all clean overlap_python

all: overlap_python

overlap_python: overlap.o fasta.o overlap_python.o
	$(CPP) $(LDFLAGS) overlap.o fasta.o overlap_python.o -o ${OUT_DIR}/coverlap.so

overlap.o: overlap.cpp fasta.h overlap.h
	$(CPP) $(CFLAGS) overlap.cpp -o overlap.o

overlap_python.o: overlap_python.cpp overlap.h
	$(CPP) $(CFLAGS) overlap_python.cpp -o overlap_python.o 

fasta.o: fasta.h fasta.cpp
	$(CPP) $(CFLAGS) fasta.cpp -o fasta.o

clean:
	rm *.o
	rm ${OUT_DIR}/coverlap.so
