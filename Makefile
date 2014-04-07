OVLP_DIR := ragout/overlap/cpp_impl/
SYNTENY_DIR := ragout/maf2synteny/cpp_impl/
PYTHON_INCLUDE := /usr/include/python2.7

ifneq ($(wildcard /usr/bin/clang),)
	CPP := clang++ -std=c++0x -stdlib=libstdc++
else
	CPP := g++ -std=c++0x
endif

.PHONY: all overlap dependencies clean maf2synteny

export PYTHON_INCLUDE
export CPP

all: overlap maf2synteny dependencies

overlap:
	make -C ${OVLP_DIR} pylib

maf2synteny:
	make -C ${SYNTENY_DIR} pylib

dependencies:
	scripts/install-deps.py

clean:
	make -C ${OVLP_DIR} pyclean
	make -C ${SYNTENY_DIR} pyclean
