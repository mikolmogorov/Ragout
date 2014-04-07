CPP := g++ -std=c++11
CFLAGS := -Wall 
LDFLAGS :=
OUT_BIN := overlap

.PHONY: all clean overlap_python debug pylib pyclean

pylib: CFLAGS += -fPIC -I${PYTHON_INCLUDE} -DPYTHON_LIB
pylib: LDFLAGS += -shared
pylib: OUT_BIN := ../coverlap.so
pylib: all

all: CFLAGS += -O3 -DNDEBUG
all: overlap_python

debug: CFLAGS += -g -D_GLIBCXX_DEBUG -D_DEBUG
debug: overlap_python

objects := ${patsubst %.cpp,%.o,${wildcard *.cpp}}
overlap_python: ${objects}
	${CPP} ${LDFLAGS} ${objects} -o ${OUT_BIN}

%.o: %.cpp *.h
	${CPP} -c ${CFLAGS} $< -o $@

pyclean:
	rm *.o
	rm ../coverlap.so

clean:
	rm *.o
	rm overlap
