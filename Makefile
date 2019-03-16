OVLP_DIR := ragout/overlap/cpp_impl
M2S_DIR := ragout/maf2synteny/cpp_impl
BIN_DIR := $(shell pwd)/bin

#setting compiler (if not set)
ifeq (${CXX},)
	IS_CLANG := $(shell which clang++ >/dev/null 2>&1; echo $$?)
	IS_GCC := $(shell which g++ >/dev/null 2>&1; echo $$?)

	ifeq (${IS_CLANG},0)
		CXX := clang++
		
	else ifeq (${IS_GCC},0)
		CXX := g++

	else
	err:
		$(error Neither gcc nor clang compilers were detected.)
	endif
endif

#adding necessary flags
CXXFLAGS += -std=c++0x
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
	CXXFLAGS += -stdlib=libc++
	LDFLAGS += -lc++
endif

export CXX
export CXXFLAGS
export LDFLAGS
export BIN_DIR

.PHONY: all overlap dependencies clean maf2synteny test

all: overlap maf2synteny

static: LDFLAGS += -static-libstdc++ -static-libgcc -static
static: all

overlap:
	make -C ${OVLP_DIR} all

maf2synteny:
	make -C ${M2S_DIR} all

dependencies:
	scripts/install-sibelia.py

test:
	scripts/run-tests.py

clean:
	make -C ${OVLP_DIR} clean
	make -C ${M2S_DIR} clean
