OVLP_DIR := ragout/overlap/cpp_impl
M2S_DIR := ragout/maf2synteny/cpp_impl
BIN_DIR := $(shell pwd)/lib

UNAME := $(shell uname -s)
IS_CLANG := $(shell which clang++ 1>&2 2>/dev/null; echo $$?)
IS_GCC := $(shell which g++ 1>&2 2>/dev/null; echo $$?)

ifeq (${IS_CLANG},0)
	CXX := clang++ -std=c++0x

	ifeq ($(UNAME),Darwin) #for macos
		CXX += -stdlib=libc++
	endif

else ifeq (${IS_GCC},0)
	CXX := g++ -std=c++0x

else
err:
	$(error Neither gcc nor clang compilers were detected.)
endif

export CXX
export BIN_DIR

.PHONY: all overlap dependencies clean maf2synteny

all: overlap maf2synteny

overlap:
	make -C ${OVLP_DIR} all

maf2synteny:
	make -C ${M2S_DIR} all

dependencies:
	scripts/install-sibelia.py

clean:
	make -C ${OVLP_DIR} clean
	make -C ${M2S_DIR} clean
