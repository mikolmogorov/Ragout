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

namespace
{
	std::vector<std::string> split(const std::string& str, const std::string& delim)
	{
		std::vector<std::string> tokens;
		size_t prev = 0, pos = 0;
		do
		{
			pos = str.find(delim, prev);
			if (pos == std::string::npos) pos = str.length();
			std::string token = str.substr(prev, pos-prev);
			if (!token.empty()) tokens.push_back(token);
			prev = pos + delim.length();
		}
		while (pos < str.length() && prev < str.length());
		return tokens;
	}
}

PermVec mafToPermutations(const std::string& mafFile, int minBlockLen)
{
	std::cerr << "\tReading maf file" << std::endl;
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

PermVec parseGff(const std::string& filename, int minBlockLen)
{
	std::cerr << "\tReading GFF file" << std::endl;

	std::unordered_map<std::string, Permutation> permBySeqId;
	//std::unordered_map<std::string, std::vector<Block>> newBlocks;

	std::ifstream fin(filename);
	if (!fin) throw std::runtime_error("Cannot open " + filename);

	std::string line;
	while(!fin.eof())
	{
		std::getline(fin, line);
		if (line.empty()) continue;
		if (line[0] == '#') continue;

		auto tokens = split(line, "\t");
		if (tokens.size() < 9) throw std::runtime_error("Error reading GFF");

		std::string seqName = tokens[0];
		int start = atoi(tokens[3].c_str());
		int end = atoi(tokens[4].c_str());
		int sign = tokens[6] == "+" ? 1 : -1;
		std::string attributes = tokens[8];
		
		int blockId = 0;
		for (auto& attr : split(attributes, ";"))
		{
			if (attr.substr(0, 2) == "id") blockId = atoi(attr.substr(3).c_str());
			//std::cout << attr << " " << blockId << std::endl;
		}
		if (blockId == 0) throw std::runtime_error("Error getting block ID");

		//std::cerr << seqName << " " << blockId << " " << start << " " << end << std::endl;
		
		permBySeqId[seqName].blocks.push_back(Block(blockId, sign, start, end));
		permBySeqId[seqName].nucLength = std::max(permBySeqId[seqName].nucLength, end);
		permBySeqId[seqName].seqName = seqName;
	}

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
