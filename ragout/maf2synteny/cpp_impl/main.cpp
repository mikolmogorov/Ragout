#ifdef PYTHON_LIB
#include <Python.h>
#include <signal.h>
#include <setjmp.h>
#endif
#include <iostream>
#include <stdexcept>

#include "breakpoint_graph.h"
#include "compress_algorithms.h"
#include "maf_tools.h"
#include "permutation.h"
#include "utility.h"

void processGraph(const PermVec& permsIn, int maxGap, PermVec& permsOut,
				  BlockGroups& groupsOut)
{
	BreakpointGraph bg(permsIn);

	DEBUG_PRINT("Started graph simplification");

	int totalPaths = 0;
	int totalBulges = 0;
	int prevPaths = 0;
	int prevBulges = 0;
	while (true)
	{
		prevPaths = compressGraph(bg, maxGap);
		DEBUG_PRINT(prevPaths << " paths compressed");
		totalPaths += prevPaths;
		if (prevPaths + prevBulges == 0) break;

		prevBulges = removeBulges(bg, maxGap);
		DEBUG_PRINT(prevBulges << " bulges removed");
		totalBulges += prevBulges;
		if (prevPaths + prevBulges == 0) break;
	}
	DEBUG_PRINT("Done: " << totalPaths << " paths compressed, "
			  	<< totalBulges << " bulges removed\n");

	bg.getPermutations(permsOut, groupsOut);
}

void compressPaths(const PermVec& permsIn, int maxGap, PermVec& permsOut,
				  BlockGroups& groupsOut)
{
	BreakpointGraph bg(permsIn);
	int paths = compressGraph(bg, maxGap);
	(void)paths;	//disable warning
	DEBUG_PRINT("Initial compression: " << paths << " paths");
	bg.getPermutations(permsOut, groupsOut);
}

struct ParamPair 
{
	int minBlock;
	int maxGap;
};

void doJob(const std::string& inputMaf, const std::string& outDir, int minBlock)
{
	std::string permsFile = outDir + "/genomes_permutations.txt";
	std::string coordsFile = outDir + "/blocks_coords.txt";
	std::string statsFile = outDir + "/coverage_report.txt";

	const int MIN_ALIGNMENT = 1;
	const int MAX_ALIGNMENT_GAP = 0;
	const auto EMPTY_GROUP = BlockGroups();
	const std::vector<ParamPair> PARAMS = {{30, 10}, {100, 100}, {500, 1000},
										  {1000, 5000}, {5000, 15000}};

	BlockGroups blockGroups;
	PermVec currentBlocks;
	PermVec mafBlocks = mafToPermutations(inputMaf, MIN_ALIGNMENT);
	compressPaths(mafBlocks, MAX_ALIGNMENT_GAP, currentBlocks, blockGroups);

	for (const ParamPair& ppair : PARAMS)
	{
		if (ppair.minBlock > minBlock) break;

		std::cerr << "Simplification with " << ppair.minBlock << " "
				  << ppair.maxGap << std::endl;
		PermVec inputBlocks = filterBySize(currentBlocks, BlockGroups(),
										   ppair.minBlock, true);
		PermVec outBlocks;
		blockGroups.clear();
		processGraph(inputBlocks, ppair.maxGap, outBlocks, blockGroups);
		if (ppair.minBlock > minBlock)
		{
			currentBlocks = mergePermutations(outBlocks, currentBlocks);
		}
		else
		{
			currentBlocks = outBlocks;
		}
	}

	PermVec outPerms = filterBySize(currentBlocks, blockGroups, 
									minBlock, true);
	renumerate(outPerms);
	outputPermutation(outPerms, permsFile);
	outputCoords(outPerms, coordsFile);
	outputStatistics(outPerms, statsFile);
}

int main(int argc, char** argv)
{
	if (argc != 4)
	{
		std::cerr << "Usage: maf2synteny <maf_file> <out_dir> <block_size>\n";
		return 1;
	}
	std::string inputMaf = argv[1];
	std::string outDir = argv[2];
	int blockSize = atoi(argv[3]);

	try
	{
		doJob(inputMaf, outDir, blockSize);
	}
	catch (std::runtime_error& e)
	{
		std::cerr << e.what() << std::endl;
		return 1;
	}

	return 0;
}

#ifdef PYTHON_LIB
static jmp_buf g_jmpEnv;
void sigintHandler(int s)
{
	longjmp(g_jmpEnv, 1);
}

static PyObject* _make_synteny(PyObject* self, PyObject* args)
{
	const char* mafFile = 0;
	const char* outDir = 0;
	int blockSize = 0;
	int terminate = -1;
	
	if (!PyArg_ParseTuple(args, "ssi", &mafFile, &outDir, &blockSize))
	{
		return Py_False;
	}

	struct sigaction pythonSig;
	sigaction(SIGINT, NULL, &pythonSig);
	signal(SIGINT, sigintHandler);

	bool result = true;
	try
	{
		terminate = setjmp(g_jmpEnv);
		if (terminate) throw std::runtime_error("SIGINT catched, exiting");
		doJob(mafFile, outDir, blockSize);
	}
	catch (std::runtime_error& e)
	{
		std::cerr << e.what() << std::endl;
		result = false;
	}
	signal(SIGINT, pythonSig.sa_handler);
	return PyBool_FromLong((long)result);
}

static PyMethodDef cmaf2syntenyMethods[] = 
{
	{"_make_synteny", _make_synteny,
	 METH_VARARGS, "Make synteny blocks from maf alignment"},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initcmaf2synteny()
{
	(void)Py_InitModule("cmaf2synteny", cmaf2syntenyMethods);
}
#endif
