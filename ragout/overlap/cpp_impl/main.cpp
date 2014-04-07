#ifdef PYTHON_LIB
#include <Python.h>
#endif

#include <iostream>

#include "overlap.h"


int main(int argc, char** argv)
{
	if (argc != 5)
	{
		std::cerr 	<< "overlap: constructs overlap graph from input contigs\n"
					<< "and outputs it in dot format\n"
					<< "Usage: overlap <fasta_in> <dot_out> <min_k> <max_k>\n";
		return 1;
	}
	return !makeOverlapGraph(argv[1], argv[2], atoi(argv[2]), atoi(argv[3]));
}

#ifdef PYTHON_LIB
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
	{"_build_overlap_graph", coverlap_build_overlap_graph,
					 METH_VARARGS, "Build overlap graph"},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initcoverlap()
{
	(void)Py_InitModule("coverlap", coverlapMethods);
}
#endif
