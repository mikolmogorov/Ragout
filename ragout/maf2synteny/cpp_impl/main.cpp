#include "breakpoint_graph.h"
#include "compress_algorithms.h"
#include "maf_tools.h"
#include "permutation.h"

#include <iostream>

void processGraph(const PermVec& permsIn, int maxGap, PermVec& permsOut,
				  BlockGroups& groupsOut)
{
	BreakpointGraph bg(permsIn);

	int totalPaths = 0;
	int totalBulges = 0;
	while (true)
	{
		int paths = compressGraph(bg, maxGap);
		int bulges = removeBulges(bg, maxGap);
		totalPaths += paths;
		totalBulges += bulges;
		if (paths + bulges == 0)
			break;
	}
	std::cerr << "Done: " << totalPaths << " paths compressed, "
			  << totalBulges << " bulges removed\n";

	bg.getPermutations(permsOut, groupsOut);
}

void compressPaths(const PermVec& permsIn, int maxGap, PermVec& permsOut,
				  BlockGroups& groupsOut)
{
	BreakpointGraph bg(permsIn);
	int paths = compressGraph(bg, maxGap);
	std::cerr << "Initial compression: " << paths << " paths\n";
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


	const int MIN_ALIGNMENT = 5;
	const int MAX_ALIGNMENT_GAP = 5;
	const float MIN_FLANK_RATE = 0.3;
	const std::vector<ParamPair> PARAMS = {{30, 30}, {100, 100}, {500, 1000},
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
										   ppair.minBlock, 0);
		PermVec outBlocks;
		blockGroups.clear();
		processGraph(inputBlocks, ppair.maxGap, outBlocks, blockGroups);
		currentBlocks = mergePermutations(outBlocks, currentBlocks, minBlock);
	}

	int flankLen = minBlock * MIN_FLANK_RATE;
	PermVec outPerms = filterBySize(currentBlocks, blockGroups, 
									minBlock, flankLen);
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
