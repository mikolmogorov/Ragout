//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

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

		if (argc == 2 && std::string(argv[1]) == "--help") return 0;
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
