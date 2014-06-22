//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <cassert>

struct Block
{
	Block(int blockId, int sign, int start, int end):
		blockId(blockId), sign(sign), start(start), end(end) {}
	int getLen() const 		{assert(end >= start); return end - start;}
	int signedId() const 	{return sign * blockId;}

	int blockId;
	int sign;
	int start;
	int end;
};

struct Permutation
{
	Permutation(int seqId=-1, const std::string& seqName="", int nucLen=-1):
		seqId(seqId), nucLength(nucLen), seqName(seqName) {}

	int seqId;
	int nucLength;
	std::vector<Block> blocks;
	std::string seqName;
};

typedef std::vector<Permutation> PermVec;
typedef std::unordered_map<int, int> BlockGroups;


inline std::unordered_map<int, Permutation*> indexBySeqId(PermVec& perms)
{
	std::unordered_map<int, Permutation*> seqId2Perm;
	for (Permutation& perm : perms)
	{
		seqId2Perm[perm.seqId] = &perm;
	}
	return seqId2Perm;
}


struct BlockPair
{
	Block* block;
	int seqId;
};

inline std::unordered_map<int, std::vector<BlockPair>> 
groupByBlockId(PermVec& perms)
{
	std::unordered_map<int, std::vector<BlockPair>> blockId2blocks;
	for (Permutation& perm : perms)
	{
		for (Block& block : perm.blocks)
		{
			blockId2blocks[block.blockId].push_back({&block, perm.seqId});
		}
	}
	return blockId2blocks;
}


void outputPermutation(const PermVec& permutations, const std::string& outFile);
void outputCoords(PermVec& permutations, const std::string& outFile);
void outputStatistics(PermVec& permutations, const std::string& outFile);
void renumerate(PermVec& permutations);
PermVec mergePermutations(PermVec& loosePerms, PermVec& finePerms);
PermVec filterBySize(PermVec& permutations, const BlockGroups& blockGroups,
					 int minBlock, bool requireAll);
