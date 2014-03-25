#include "BreakpointGraph.h"
#include "Utility.h"
#include <iostream>

BreakpointGraph::BreakpointGraph(const std::vector<Permutation>& permutations)
{
	for (auto &perm : permutations)
	{
		//black edges
		for (auto &block : perm.blocks)
		{
			std::vector<Edge*> edges;
			this->getBlackEdges(block.blockId, -block.blockId, edges);
			if (edges.empty())
				this->addEdge(block.blockId, -block.blockId, Edge::BLACK);
		}

		//chromosome ends
		Edge* headEdge = this->addEdge(INFINUM, perm.blocks.front().blockId,
									   perm.seqId);
		headEdge->rightPos = perm.blocks.front().start;
		Edge* tailEdge = this->addEdge(-perm.blocks.back().blockId, INFINUM,
									   perm.seqId);
		tailEdge->leftPos = perm.blocks.back().end;
		_origins.push_back(headEdge);

		Edge* prevEdge = headEdge;
		Edge* curEdge = headEdge;

		//adjacencies
		for (auto blockPair : make_adjacent_range(perm.blocks))
		{
			int leftPos = blockPair.first.end;
			int rightPos = blockPair.second.start;

			if (rightPos < leftPos)
				std::cerr << "WRANING: overlapping blocks\n"
						  << blockPair.first.start << " " << blockPair.first.end
						  << " " << blockPair.second.start << " "
						  << blockPair.second.end << " | " << perm.seqId << "\n";

			curEdge = this->addEdge(-blockPair.first.blockId,
									blockPair.second.blockId, perm.seqId);
			curEdge->leftPos = leftPos;
			curEdge->rightPos = rightPos;

			prevEdge->nextEdge = curEdge;
			curEdge->prevEdge = prevEdge;
			prevEdge = curEdge;
		}

		curEdge->nextEdge = tailEdge;
		tailEdge->prevEdge = curEdge;
	}
}

void BreakpointGraph::getFragmentedBlocks(std::vector<std::vector<int>>& groups)
{
	//TODO
}

void BreakpointGraph::getPermutations(std::vector<Permutation>& permutations)
{
	int nextEdgeId = 0;
	std::unordered_map<Edge*, int> edgeIds;
	auto getEdgeId = [&nextEdgeId, &edgeIds] (Edge* e) 
	{
		if (!edgeIds.count(e)) edgeIds[e] = nextEdgeId++;
		return edgeIds[e];
	};
	///
	
	for (auto startEdge : _origins)
	{
		permutations.push_back(Permutation());
		Edge* prevEdge = startEdge;
		Edge* curEdge = startEdge->nextEdge;

		while (curEdge)
		{
			std::vector<Edge*> blackEdges;
			this->getBlackEdges(prevEdge->rightNode, curEdge->leftNode, 
								blackEdges);
			assert(blackEdges.size() == 1);

			int blockId = getEdgeId(blackEdges.back());
			int sign = (blackEdges.back()->rightNode == curEdge->leftNode) 
							? 1 : -1;
			int start = prevEdge->rightPos;
			int end = curEdge->leftPos;
			assert(end >= start);

			permutations.back().blocks.push_back(Block(blockId, sign, 
													   start, end));
			permutations.back().seqId = startEdge->seqId;
			prevEdge = curEdge;
			curEdge = curEdge->nextEdge;
		}
	}
}
