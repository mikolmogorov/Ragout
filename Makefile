PY=python
OVLP_DIR="ragout/overlap/cpp_impl/"
export PYTHON_INCLUDE="/usr/include/python2.7"

all: overlap dependencies

overlap:
	cd ${OVLP_DIR}; make

dependencies:
	${PY} scripts/install-deps.py

clean:
	cd ${OVLP_DIR}; make clean
