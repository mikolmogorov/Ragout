OVLP_DIR := ragout/overlap/cpp_impl/
M2S_DIR := ragout/maf2synteny/cpp_impl/
BIN_DIR := $(shell pwd)/lib/

UNAME := $(shell uname -s)
ifneq ($(wildcard /usr/bin/clang),)
	CPP := clang++ -std=c++11

	ifeq ($(UNAME),Darwin) #for macos
		CPP += -stdlib=libc++
	endif
else
	CPP := g++ -std=c++11
endif

export CPP
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