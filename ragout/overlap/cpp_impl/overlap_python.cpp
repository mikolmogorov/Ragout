#include <Python.h>
#include <iostream>
#include "overlap.h"


static PyObject*
coverlap_build_overlap_graph(PyObject* self, PyObject* args)
{
	const char* fileIn = 0;
	const char* fileOut = 0;
	int minOverlap = 0;
	int maxOverlap = 0;

	if (!PyArg_ParseTuple(args, "ssii", &fileIn, &fileOut,
						  &minOverlap, &maxOverlap))
		return Py_False;
	bool ret = makeOverlapGraph(fileIn, fileOut, minOverlap, maxOverlap);
	return PyBool_FromLong((long)ret);
}

static PyMethodDef coverlapMethods[] = 
{
	{"build_overlap_graph", coverlap_build_overlap_graph,
					METH_VARARGS, "Build overlap graph"},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initcoverlap()
{
	(void)Py_InitModule("coverlap", coverlapMethods);
}
