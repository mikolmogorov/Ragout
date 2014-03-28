#include "breakpoint_graph.h"
#include "compress_algorithms.h"
#include "maf_tools.h"
#include "permutation.h"

#include <iostream>

void simplifyGraph(BreakpointGraph& bg, int maxGap)
{
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

	std::string permsFile = outDir + "/genomes_permutations.txt";

	std::vector<Permutation> perms = mafToPermutations(inputMaf, 100);
	BreakpointGraph bg(perms);
	simplifyGraph(bg, 1000);
	PermVec outPerms = bg.getPermutations();

	outputPermutation(outPerms, permsFile);
	return 0;
}
