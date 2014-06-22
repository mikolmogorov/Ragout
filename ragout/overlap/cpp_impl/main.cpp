//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#ifdef PYTHON_LIB
#include <Python.h>
#include <signal.h>
#include <setjmp.h>
#endif

#include <iostream>

#include "build_graph.h"


int main(int argc, char** argv)
{
	if (argc != 5)
	{
		std::cerr 	<< "overlap: constructs overlap graph from input contigs\n"
					<< "and outputs it in dot format\n"
					<< "Usage: overlap <fasta_in> <dot_out> <min_k> <max_k>\n";
		return 1;
	}
	return !makeOverlapGraph(argv[1], argv[2], atoi(argv[3]),
							 atoi(argv[4]), true);
}

#ifdef PYTHON_LIB

static jmp_buf g_jumpEnv;

void sigintHandler(int sig)
{
	longjmp(g_jumpEnv, 1);
	//_exit(1);
}

static PyObject*
coverlap_build_overlap_graph(PyObject* self, PyObject* args)
{
	const char* fileIn = 0;
	const char* fileOut = 0;
	int minOverlap = 0;
	int maxOverlap = 0;
	int terminate = -1;
	bool filterKmer = true;

	if (!PyArg_ParseTuple(args, "ssiib", &fileIn, &fileOut,
						  &minOverlap, &maxOverlap, &filterKmer))
	{
		return Py_False;
	}

	struct sigaction pythonSig;
	sigaction(SIGINT, NULL, &pythonSig);
	signal(SIGINT, sigintHandler);

	bool result;
	terminate = setjmp(g_jumpEnv);
	if (!terminate)
	{
		result = makeOverlapGraph(fileIn, fileOut, minOverlap,
								  maxOverlap, filterKmer);
	}
	else
	{
		std::cerr << "SIGINT catched, exiting\n";
		result = false;
	}
	signal(SIGINT, pythonSig.sa_handler);
	return PyBool_FromLong((long)result);
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
