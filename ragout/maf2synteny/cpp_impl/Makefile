CFLAGS += -Wall
OUT_BIN := ../maf2synteny

.PHONY: all clean debug pylib pyclean

pylib: CFLAGS += -fPIC -I${PYTHON_INCLUDE} -DPYTHON_LIB
pylib: LDFLAGS += -shared
pylib: OUT_BIN := ../cmaf2synteny.so
pylib: all

all: CFLAGS += -O2 -DNDEBUG
all: maf2synteny

debug: CFLAGS += -g -D_GLIBCXX_DEBUG -D_DEBUG
debug: maf2synteny

objects := ${patsubst %.cpp,%.o,${wildcard *.cpp}}
maf2synteny: ${objects}
	${CPP} ${LDFLAGS} ${objects} -o ${OUT_BIN}

%.o: %.cpp *.h
	${CPP} -c ${CFLAGS} $< -o $@

pyclean:
	rm *.o
	rm ../cmaf2synteny.so

clean:
	rm *.o
	rm ../maf2synteny
