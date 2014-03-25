OVLP_DIR := ragout/overlap/cpp_impl/
PYTHON_INCLUDE := /usr/include/python2.7

ifneq ($(wildcard /usr/bin/clang),)
	CPP := clang++ -std=c++0x -stdlib=libstdc++
else
	CPP := g++ -std=c++0x
endif

.PHONY: all overlap dependencies clean

export PYTHON_INCLUDE
export CPP

all: overlap dependencies

overlap:
	make -C ${OVLP_DIR}

dependencies:
	scripts/install-deps.py

clean:
	make clean -C ${OVLP_DIR}
