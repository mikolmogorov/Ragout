#include "permutation.h"
#include <stdexcept>
#include <fstream>

void outputPermutation(const PermVec& permutations, const std::string outFile)
{
	std::ofstream fout(outFile);
	if (!fout)
		throw std::runtime_error("Can't open " + outFile);
	for (const Permutation& perm : permutations)
	{
		fout << ">" << perm.seqName << std::endl;
		for (const Block& block : perm.blocks)
		{
			char sign = (block.sign > 0) ? '+' : '-';
			fout << sign << block.blockId << " ";
		}
		fout << "$\n";
	}
}
