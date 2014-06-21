CXXFLAGS += -Wall
OUT_BIN := ${BIN_DIR}/ragout-overlap

.PHONY: all clean debug

all: CXXFLAGS += -O2 -DNDEBUG
all: maf2synteny

debug: CXXFLAGS += -g -D_GLIBCXX_DEBUG -D_DEBUG
debug: maf2synteny

objects := ${patsubst %.cpp,%.o,${wildcard *.cpp}}
maf2synteny: ${objects}
	${CXX} ${LDFLAGS} ${objects} -o ${OUT_BIN}

%.o: %.cpp *.h
	${CXX} -c ${CXXFLAGS} $< -o $@

clean:
	rm *.o
	rm ${OUT_BIN}
