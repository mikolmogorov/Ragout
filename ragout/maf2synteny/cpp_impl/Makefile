CFLAGS += -Wall
OUT_BIN := ${BIN_DIR}/ragout-maf2synteny

.PHONY: all clean debug

all: CFLAGS += -O2 -DNDEBUG
all: maf2synteny

debug: CFLAGS += -g -D_GLIBCXX_DEBUG -D_DEBUG
debug: maf2synteny

objects := ${patsubst %.cpp,%.o,${wildcard *.cpp}}
maf2synteny: ${objects}
	${CPP} ${LDFLAGS} ${objects} -o ${OUT_BIN}

%.o: %.cpp *.h
	${CPP} -c ${CFLAGS} $< -o $@

clean:
	rm *.o
	rm ${OUT_BIN}
