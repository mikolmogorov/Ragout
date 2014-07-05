//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#ifdef PYTHON_LIB
#include <Python.h>
#include <signal.h>
#include <setjmp.h>
#endif
#include <iostream>
#include <stdexcept>
#include <fstream>

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
	std::cerr << "\tStarted initial compression\n";
	BreakpointGraph bg(permsIn);
	int paths = compressGraph(bg, maxGap);
	(void)paths;	//disable warning
	DEBUG_PRINT("Initial compression: " << paths << " paths");
	bg.getPermutations(permsOut, groupsOut);
}

void outputBlocks(const std::string& outDir, PermVec& blocks,
				  BlockGroups& groups, int blockSize)
{
	makeDirectory(outDir);
	std::string permsFile = outDir + "/genomes_permutations.txt";
	std::string coordsFile = outDir + "/blocks_coords.txt";
	std::string statsFile = outDir + "/coverage_report.txt";

	PermVec outPerms = filterBySize(blocks, groups, blockSize, true);
	renumerate(outPerms);
	outputPermutation(outPerms, permsFile);
	outputCoords(outPerms, coordsFile);
	outputStatistics(outPerms, statsFile);
}

struct ParamPair 
{
	int minBlock;
	int maxGap;
};


std::vector<ParamPair> parseSimplParamsFile(const std::string& filename)
{
	std::vector<ParamPair> out;
	std::ifstream fin(filename);
	if (!fin) throw std::runtime_error("Could not open " + filename);

	std::string buffer;
	while(!fin.eof())
	{
		std::getline(fin, buffer);
		if (buffer.empty()) break;

		size_t sep = buffer.find_first_of(" \t");
		if (sep == std::string::npos)
		{
			throw std::runtime_error("Error parsing " + filename);
		}
		
		int k = atoi(buffer.substr(0, sep).c_str());
		int d = atoi(buffer.substr(sep + 1).c_str());
		out.push_back({k, d});
	}
	return out;
}

void doJob(const std::string& inputMaf, const std::string& outDir, 
		   std::vector<ParamPair> simplParams, std::vector<int> minBlockSizes)
{

	const int MIN_ALIGNMENT = 1;
	const int MAX_ALIGNMENT_GAP = 0;
	const auto EMPTY_GROUP = BlockGroups();

	BlockGroups blockGroups = EMPTY_GROUP;
	PermVec currentBlocks;
	//sort blocks in reverse order (will use it as stack)
	std::sort(minBlockSizes.begin(), minBlockSizes.end(), std::greater<int>());
	makeDirectory(outDir);

	//read maf alignment and join adjacent columns
	std::cerr << "\tParsing MAF file\n";
	PermVec mafBlocks = mafToPermutations(inputMaf, MIN_ALIGNMENT);
	compressPaths(mafBlocks, MAX_ALIGNMENT_GAP, currentBlocks, blockGroups);

	//iterative simplification
	for (const ParamPair& ppair : simplParams)
	{
		if (minBlockSizes.empty()) break;
		//output blocks of certain size
		while (!minBlockSizes.empty() && minBlockSizes.back() < ppair.minBlock)
		{
			std::string blockDir = outDir + "/" + 
								   std::to_string(minBlockSizes.back());
			outputBlocks(blockDir, currentBlocks, blockGroups,
						 minBlockSizes.back());
			minBlockSizes.pop_back();
		}

		std::cerr << "\tSimplification with " << ppair.minBlock << " "
				  << ppair.maxGap << std::endl;
		PermVec inputBlocks = filterBySize(currentBlocks, BlockGroups(),
										   ppair.minBlock, true);
		PermVec outBlocks;
		blockGroups.clear();
		processGraph(inputBlocks, ppair.maxGap, outBlocks, blockGroups);
		currentBlocks = outBlocks;
	}

	//if any left
	for (int minBlock : minBlockSizes)
	{
		std::string blockDir = outDir + "/" + std::to_string(minBlock);
		outputBlocks(blockDir, currentBlocks, blockGroups, minBlock);
	}
}

int main(int argc, char** argv)
{
	if (argc < 4)
	{
		std::cerr << "Usage: maf2synteny maf_file out_dir simpl_params "
				  << "block_size_1 [block_size_2 ...]\n";
		return 1;
	}
	std::string inputMaf = argv[1];
	std::string outDir = argv[2];
	std::string paramsFile = argv[3];

	std::vector<int> sizes;
	for (int i = 4; i < argc; ++i)
	{
		sizes.push_back(atoi(argv[i]));
	}

	try
	{
		std::vector<ParamPair> simplParams = parseSimplParamsFile(paramsFile);
		doJob(inputMaf, outDir, simplParams, sizes);
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
	std::vector<int> blockSizes;
	int terminate = -1;

	PyObject* listObj = nullptr;
	if (!PyArg_ParseTuple(args, "ssO", &mafFile, &outDir, &listObj))
	{
		return Py_False;
	}
	int nBlocks = PyList_Size(listObj);
	for (int i = 0; i < nBlocks; ++i)
	{
		blockSizes.push_back(PyInt_AsLong(PyList_GetItem(listObj, i)));
	}

	struct sigaction pythonSig;
	sigaction(SIGINT, NULL, &pythonSig);
	signal(SIGINT, sigintHandler);

	bool result = true;
	try
	{
		terminate = setjmp(g_jmpEnv);
		if (terminate) throw std::runtime_error("SIGINT catched, exiting");
		doJob(mafFile, outDir, blockSizes);
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
