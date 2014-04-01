#CPP := g++ -std=c++11 -g -D_GLIBCXX_DEBUG
CPP := g++ -std=c++11 -O3
CFLAGS := -c -Wall
LDFLAGS := 

.PHONY: all clean

all: maf2synteny

objects := ${patsubst %.cpp,%.o,${wildcard *.cpp}}
maf2synteny: ${objects}
	${CPP} ${LDFLAGS} ${objects} -o $@

%.o: %.cpp *.h
	${CPP} ${CFLAGS} $< -o $@

clean:
	rm *.o
	rm maf2synteny
