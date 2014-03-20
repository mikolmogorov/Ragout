PY=python
LIB_DIR=lib

all: overlap dependencies

overlap:
	cd ${LIB_DIR}/overlap/src; make

dependencies:
	${PY} bin/install-deps.py

clean:
	cd ${LIB_DIR}/overlap/src; make clean
