#include "permutation.h"
#include <stdexcept>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

static const std::string SEPARATOR(80, '-');

void outputPermutation(const PermVec& permutations, const std::string& outFile)
{
	std::ofstream fout(outFile);
	if (!fout) throw std::runtime_error("Can't open " + outFile);
	for (const Permutation& perm : permutations)
	{
		fout << ">" << perm.seqName << std::endl;
		for (const Block& block : perm.blocks)
		{
			char strand = (block.sign > 0) ? '+' : '-';
			fout << strand << block.blockId << " ";
		}
		fout << "$\n";
	}
}

void outputCoords(const PermVec& permutations, const std::string& outFile)
{
	std::ofstream fout(outFile);
	if (!fout) throw std::runtime_error("Can't open " + outFile);

	std::unordered_map<int, std::vector<const Block*>> byBlock;
	std::unordered_map<const Block*, int> seqIds;
	for (const Permutation& perm : permutations)
	{
		for (const Block& block : perm.blocks)
		{
			byBlock[block.blockId].push_back(&block);
			seqIds[&block] = perm.seqId;
		}
	}
	
	fout << "Seq_id\tSize\tDescription\n";
	for (const Permutation& perm : permutations)
	{
		fout << perm.seqId << "\t" << perm.nucLength << "\t"
			 << perm.seqName << std::endl;
	}
	fout << SEPARATOR << std::endl;

	for (auto &itBlocks : byBlock)
	{
		fout << "Block #" << itBlocks.first << "\nSeq_id\tStrand\tStart\t"
			 << "End\tLength\n";

		for (const Block* block : itBlocks.second)
		{
			char strand = (block->sign > 0) ? '+' : '-';
			fout << seqIds[block] << "\t" << strand << "\t" << block->start
				 << "\t" << block->end << "\t" << block->getLen() << std::endl;
		}
		fout << SEPARATOR << std::endl;
	}
}

void outputStatistics(PermVec& permutations, const std::string& outFile)
{
	std::ofstream fout(outFile);
	if (!fout) throw std::runtime_error("Can't open " + outFile);

	std::unordered_map<int, std::vector<const Block*>> byBlock;
	std::unordered_map<int, int> multiplicity;
	std::unordered_map<std::string, float> covered;

	for (const Permutation& perm : permutations)
	{
		for (const Block& block : perm.blocks)
		{
			byBlock[block.blockId].push_back(&block);
			assert(!perm.seqName.empty());
			covered[perm.seqName] += block.getLen();
		}
		covered[perm.seqName] /= perm.nucLength;
	}
	
	fout << "Seq_id\tSize\tDescription\n";
	for (const Permutation& perm : permutations)
	{
		fout << perm.seqId << "\t" << perm.nucLength << "\t"
			 << perm.seqName << std::endl;
	}
	fout << SEPARATOR << std::endl;

	for (auto &blockPair : byBlock)
	{
		++multiplicity[blockPair.second.size()];
	}
	for (auto &mulPair : multiplicity)
	{
		fout << mulPair.first << "\t" << mulPair.second << std::endl;
	}
	fout << SEPARATOR << std::endl;

	for (auto &covPair : covered)
	{
		fout << covPair.first << "\t" << covPair.second * 100 << std::endl;
	}
}

void renumerate(PermVec& permutations)
{
	int nextId = 1;
	std::unordered_map<int, int> newIds;
	auto newId = [&nextId, &newIds] (int oldId)
	{
		if (!newIds.count(oldId))
			newIds[oldId] = nextId++;
		return newIds[oldId];
	};

	for (Permutation& perm : permutations)
	{
		for (Block& block : perm.blocks)
		{
			block.blockId = newId(block.blockId);
		}
	}
}

//the function merges two permutations in different scales.
//simplifiedPerms is expoected to be in finer, 
//and initialPerms - in lower;
PermVec mergePermutations(const PermVec& simplifiedPerms,
						  const PermVec& initialPerms)
{
	PermVec permutations = simplifiedPerms;
	std::unordered_map<int, Permutation*> bySeqId;
	for (Permutation& perm : permutations)
	{
		bySeqId[perm.seqId] = &perm;
	}

	std::unordered_map<int, std::vector<int>> blockStarts;
	std::unordered_map<int, std::vector<int>> blockEnds;
	for (const Permutation& perm : simplifiedPerms)
	{
		for (const Block& block : perm.blocks)
		{
			blockStarts[perm.seqId].push_back(block.start);
			blockEnds[perm.seqId].push_back(block.end);
		}
	}

	int nextId = 1;
	std::unordered_map<int, std::vector<const Block*>> byBlock;
	std::unordered_map<const Block*, int> seqIds;
	for (const Permutation& perm : initialPerms)
	{
		for (const Block& block : perm.blocks)
		{
			byBlock[block.blockId].push_back(&block);
			seqIds[&block] = perm.seqId;
			nextId = std::max(nextId, block.blockId);
		}
	}
	++nextId;

	//here we check if block from finer scale do not intersect with
	//others from loose scale
	for (auto &blockPair : byBlock)
	{
		std::vector<const Block*> toInsert;
		for (const Block* block : blockPair.second)
		{
			auto &endVec = blockEnds[seqIds[block]];
			auto leftIns = std::upper_bound(endVec.begin(), endVec.end(), 
										   block->start);
			auto &startVec = blockStarts[seqIds[block]];
			auto rightIns = std::upper_bound(startVec.begin(), startVec.end(), 
										    block->end);

			if (leftIns == rightIns) toInsert.push_back(block);
		}

		if (!toInsert.empty())
		{
			for (const Block* block : toInsert)
			{
				int seqId = seqIds[block];
				bySeqId[seqId]->blocks.push_back(*block);
				bySeqId[seqId]->blocks.back().blockId = nextId;
			}
			++nextId;
		}
	}

	auto cmp = [](const Block& a, const Block& b) {return a.start < b.start;};
	for (Permutation& perm : permutations)
	{
		std::sort(perm.blocks.begin(), perm.blocks.end(), cmp);
	}

	return permutations;
}

PermVec filterBySize(const PermVec& permutations, 
					 const BlockGroups& blockGroups, int minBlock, int minFlank)
{
	PermVec outPerms;

	std::unordered_map<int, std::unordered_map<int, int>> groupLen;
	for (const Permutation& perm : permutations)
	{
		for (const Block& block : perm.blocks)
		{
			assert(block.blockId);
			if (blockGroups.count(block.blockId))
			{
				int groupId = blockGroups.at(block.blockId);
				groupLen[perm.seqId][groupId] += block.getLen();
			}
		}
	}

	std::unordered_set<int> shouldOutput;
	for (const Permutation& perm : permutations)
	{
		for (const Block& block : perm.blocks)
		{
			if (block.getLen() >= minBlock)
			{
				shouldOutput.insert(block.blockId);
			}
			else
			{
				auto groupId = blockGroups.find(block.blockId);
				if (groupId != blockGroups.end() &&
						groupLen[perm.seqId][groupId->second] >= minBlock &&
						block.getLen() >= minFlank)
					shouldOutput.insert(block.blockId);
			}
		}
	}

	for (const Permutation& perm : permutations)
	{
		outPerms.push_back(Permutation());
		outPerms.back().seqId = perm.seqId;
		outPerms.back().nucLength = perm.nucLength;
		outPerms.back().seqName = perm.seqName;
		for (const Block& block : perm.blocks)
		{
			if (shouldOutput.count(block.blockId))
			{
				outPerms.back().blocks.push_back(block);
			}
		}
		if (outPerms.back().blocks.empty()) outPerms.pop_back();
	}

	return outPerms;
}


