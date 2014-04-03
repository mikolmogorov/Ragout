CPP := g++ -std=c++11
CFLAGS := -Wall
LDFLAGS := 

.PHONY: all clean debug

all: CFLAGS += -O3
all: maf2synteny

debug: CFLAGS += -g -D_GLIBCXX_DEBUG -D_DEBUG
debug: maf2synteny

objects := ${patsubst %.cpp,%.o,${wildcard *.cpp}}
maf2synteny: ${objects}
	${CPP} ${LDFLAGS} ${objects} -o $@

%.o: %.cpp *.h
	${CPP} -c ${CFLAGS} $< -o $@

clean:
	rm *.o
	rm maf2synteny
