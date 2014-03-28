#pragma once

#include <vector>
#include <string>
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

void outputPermutation(const PermVec& permutations, const std::string outFile);
