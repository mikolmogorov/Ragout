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
	int seqId;
	int nucLength;
	std::vector<Block> blocks;
	std::string seqName;
};

typedef std::vector<Permutation> PermVec;
typedef std::unordered_map<int, int> BlockGroups;

void outputPermutation(const PermVec& permutations, const std::string& outFile);
void outputCoords(const PermVec& permutations, const std::string& outFile);
void outputStatistics(PermVec& permutations, const std::string& outFile);
void renumerate(PermVec& permutations);
PermVec mergePermutations(const PermVec& simplifiedPerms,
						  const PermVec& initialPerms);
PermVec filterBySize(const PermVec& permutations, 
					 const BlockGroups& blockGroups, int minBlock, int minFlank);
