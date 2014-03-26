#include "MafTools.h"
#include "Utility.h"
#include <fstream>
#include <stdexcept>
#include <unordered_map>

std::vector<Permutation> mafToPermutations(const std::string& mafFile, 
										   int minBlockLen)
{
	const int MAX_GAP_RATE = 0.3;

	int blockId = 0;
	std::unordered_map<std::string, Permutation> permBySeqId;
	std::unordered_map<std::string, Permutation> newBlocks;
	bool newLcb = false;

	auto updatePerms = [&permBySeqId, &newBlocks] ()
	{
		if (newBlocks.size() > 1)
			for (auto p : newBlocks)
				std::copy(p.second.blocks.begin(), 
						  p.second.blocks.end(),
						  std::back_inserter(permBySeqId[p.first].blocks));
		newBlocks.clear();
	};

	std::ifstream fin(mafFile);
	if (!fin)
		throw std::runtime_error("Cannot open " + mafFile);
	std::string line;
	while(!fin.eof())
	{
		std::getline(fin, line);
		if (line.empty())
			continue;

		if (line[0] == 'a')
		{
			newLcb = true;
			updatePerms();
		}
		else if (line[0] == 's')
		{
			std::vector<std::string> tokens = split(line, '\t');
		}
	}

	std::vector<Permutation> permutations;
}
