//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#include "maf_tools.h"
#include "utility.h"

#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <sstream>
#include <iostream>

PermVec mafToPermutations(const std::string& mafFile, int minBlockLen)
{
	DEBUG_PRINT("Reading maf file");
	const float MAX_GAP_RATE = 0.3;

	int blockId = 1;
	std::unordered_map<std::string, Permutation> permBySeqId;
	std::unordered_map<std::string, std::vector<Block>> newBlocks;

	auto updatePerms = [&permBySeqId, &newBlocks, &blockId] ()
	{
		if (newBlocks.size() > 1)
		{
			for (auto p : newBlocks)
				std::copy(p.second.begin(), p.second.end(),
						  std::back_inserter(permBySeqId[p.first].blocks));
			++blockId;
		}
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
			updatePerms();
		}
		else if (line[0] == 's')
		{
			std::stringstream ss(line.substr(2));
			std::string seqName;
			std::string strand;
			std::string seq;
			int start = -1;
			int srcLen = -1;
			int ungappedLen = -1;
			ss >> seqName >> start >> ungappedLen >> strand >> srcLen >> seq;

			//adhoc fix for progressiveCactus
			size_t dotPos = seqName.find('.');
			std::string chrName = (dotPos != seqName.npos) ? 
								  seqName.substr(dotPos + 1) : seqName;
			if (chrName.substr(0, 2) == "gi" && seqName.back() != '|')
			{
				seqName += "|";
			}
			//
			
			int absoluteStart = (strand == "+") ? start : 
								srcLen - (start + ungappedLen);
			float gapRate = float(seq.length() - ungappedLen) / seq.length();
			if (ungappedLen >= minBlockLen && gapRate < MAX_GAP_RATE)
			{
				int sign = (strand == "+") ? 1 : -1;
				newBlocks[seqName].push_back(Block(blockId, sign, absoluteStart,
												 absoluteStart + ungappedLen));
			}

			permBySeqId[seqName].nucLength = srcLen;
			permBySeqId[seqName].seqName = seqName;
		}
	}
	updatePerms();

	std::vector<Permutation> permutations;
	int seqId = 1;
	auto cmp = [](const Block& a, const Block& b) {return a.start < b.start;};
	for (auto p : permBySeqId)
	{
		if (!p.second.blocks.empty())
		{
			permutations.push_back(std::move(p.second));
			permutations.back().seqId = seqId++;
			std::sort(permutations.back().blocks.begin(), 
					  permutations.back().blocks.end(), cmp);
		}
	}
	DEBUG_PRINT("Finished reading");
	return permutations;
}
