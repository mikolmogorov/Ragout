//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#include "permutation.h"
#include "utility.h"
#include <stdexcept>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>

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

void outputCoords(PermVec& permutations, const std::string& outFile)
{
	std::ofstream fout(outFile);
	if (!fout) throw std::runtime_error("Can't open " + outFile);

	fout << "Seq_id\tSize\tDescription\n";
	for (Permutation& perm : permutations)
	{
		fout << perm.seqId << "\t" << perm.nucLength << "\t"
			 << perm.seqName << std::endl;
	}
	fout << SEPARATOR << std::endl;

	auto blocksIndex = groupByBlockId(permutations);
	for (auto &itBlocks : blocksIndex)
	{
		fout << "Block #" << itBlocks.first << "\nSeq_id\tStrand\tStart\t"
			 << "End\tLength\n";

		for (BlockPair& bp : itBlocks.second)
		{
			char strand = (bp.block->sign > 0) ? '+' : '-';
			int convStart = bp.block->start;
			int convEnd = bp.block->end;
			if (strand == '-') std::swap(convStart, convEnd);

			fout << bp.seqId << "\t" << strand << "\t"
				 << convStart << "\t" << convEnd
				 << "\t" << bp.block->getLen() << std::endl;
		}
		fout << SEPARATOR << std::endl;
	}
}

void outputStatistics(PermVec& permutations, const std::string& outFile)
{
	std::ofstream fout(outFile);
	if (!fout) throw std::runtime_error("Can't open " + outFile);

	std::unordered_map<int, int> multiplicity;
	std::unordered_map<std::string, float> covered;

	for (const Permutation& perm : permutations)
	{
		for (const Block& block : perm.blocks)
		{
			covered[perm.seqName] += block.getLen();
		}
		covered[perm.seqName] /= perm.nucLength;
	}
	
	//fout << "Seq_id\tSize\tDescription\n";
	//for (const Permutation& perm : permutations)
	//{
	//	fout << perm.seqId << "\t" << perm.nucLength << "\t"
	//		 << perm.seqName << std::endl;
	//}
	//fout << SEPARATOR << std::endl;

	auto blocksIndex = groupByBlockId(permutations);
	for (auto &blockPair : blocksIndex)
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
	std::unordered_map<int, int> newSigns;
	auto newIdAndSign = [&nextId, &newIds, &newSigns] (Block& block)
	{
		if (!newIds.count(block.blockId))
		{
			newIds[block.blockId] = nextId++;
			newSigns[block.blockId] = block.sign;
		}
		return std::make_pair(newIds[block.blockId], newSigns[block.blockId]);
	};

	for (Permutation& perm : permutations)
	{
		for (Block& block : perm.blocks)
		{
			auto idAndSign = newIdAndSign(block);
			block.blockId = idAndSign.first;
			block.sign *= idAndSign.second;
		}
	}
}

//the function merges two permutations in different scales.
PermVec mergePermutations(PermVec& loosePerms, PermVec& finePerms)
{
	std::unordered_map<int, std::vector<int>> blockStarts;
	std::unordered_map<int, std::vector<int>> blockEnds;
	int nextId = 0;
	for (Permutation& perm : loosePerms)
	{
		for (Block& block : perm.blocks)
		{
			blockStarts[perm.seqId].push_back(block.start);
			blockEnds[perm.seqId].push_back(block.end);
			nextId = std::max(nextId, block.blockId);
		}
	}
	++nextId;

	//here we check if block from finer scale do not intersect with
	//others from loose scale
	auto fineIndex = groupByBlockId(finePerms);
	std::vector<int> blocksToInsert;
	for (auto &indexPair : fineIndex)
	{
		bool inserting = true;
		for (BlockPair& bp : indexPair.second)
		{
			auto &endVec = blockEnds[bp.seqId];
			int leftIns = std::upper_bound(endVec.begin(), endVec.end(),
										   bp.block->start) - endVec.begin();
			auto &startVec = blockStarts[bp.seqId];
			int rightIns = std::upper_bound(startVec.begin(), startVec.end(),
											bp.block->end) - startVec.begin();

			if (leftIns != rightIns) 
			{
				inserting = false;
				break;
			}
		}

		if (inserting) blocksToInsert.push_back(indexPair.first);
	}

	std::unordered_map<int, std::vector<Block>> outBlocks;
	auto fineBySeqId = indexBySeqId(finePerms);
	for (Permutation& perm : loosePerms)
	{
		outBlocks[perm.seqId] = perm.blocks;
	}
	for (int block : blocksToInsert)
	{
		for (BlockPair& bp : fineIndex[block])
		{
			outBlocks[bp.seqId].push_back(*bp.block);
			outBlocks[bp.seqId].back().blockId = nextId;
		}
		++nextId;
	}

	PermVec outPerms;
	auto cmp = [](const Block& a, const Block& b) {return a.start < b.start;};
	for (auto &blockPair : outBlocks)
	{
		std::sort(blockPair.second.begin(), blockPair.second.end(), cmp);

		Permutation* p = fineBySeqId[blockPair.first];
		outPerms.push_back(Permutation(blockPair.first, p->seqName, 
									   p->nucLength));
		outPerms.back().blocks = std::move(blockPair.second);
	}

	return outPerms;
}

PermVec filterBySize(PermVec& permutations, const BlockGroups& blockGroups,
					 int minBlock, bool requireAll)
{
	//const float MIN_FLANK_RATE = 0.3;
	//const int minFlank = minBlock * MIN_FLANK_RATE;

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
	auto blocksIndex = groupByBlockId(permutations);
	//int acceptedGroup = 0;
	for (auto &itBlocks : blocksIndex)
	{
		size_t numAccepted = 0;
		for (BlockPair& bp : itBlocks.second)
		{
			if (bp.block->getLen() >= minBlock)
			{
				++numAccepted;
			}
			//else
			//{
			//	auto groupId = blockGroups.find(bp.block->blockId);
			//	if (groupId != blockGroups.end() &&
			//		groupLen[bp.seqId][groupId->second] >= minBlock &&
			//		bp.block->getLen() >= minFlank)
			//	{
			//		++numAccepted;
			//		++acceptedGroup;
			//	}
			//}
		}
		if (numAccepted == itBlocks.second.size() ||
			(numAccepted && !requireAll))
		{
			shouldOutput.insert(itBlocks.first);
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
	//DEBUG_PRINT("Accepted by group: " << acceptedGroup);

	return outPerms;
}
