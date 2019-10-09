//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#include "breakpoint_graph.h"
#include "utility.h"
#include "compress_algorithms.h"
#include "disjoint_set.h"

#include <unordered_set>
#include <iostream>
#include <set>

BreakpointGraph::~BreakpointGraph()
{
	DEBUG_PRINT("Graph destructor");

	std::unordered_set<Edge*> toDelete;
	for (int node : this->iterNodes())
	{
		//this->removeNode(node);
		for (Edge* e : _nodes[node].edges) toDelete.insert(e);
	}
	for (Edge* e : toDelete) delete e;

	DEBUG_PRINT("Graph destructor - finished");
}

BreakpointGraph::BreakpointGraph(const std::vector<Permutation>& permutations)
{
	DEBUG_PRINT("Building breakpoint graph");
	int maxBlockOverlap = 0;
	for (auto &perm : permutations)
	{
		assert(!perm.blocks.empty());
		_seqNames[perm.seqId] = perm.seqName;
		_seqLength[perm.seqId] = perm.nucLength;

		//black edges
		for (auto &block : perm.blocks)
		{
			if (this->getBlackEdges(block.blockId, -block.blockId).empty())
			{
				this->addEdge(block.blockId, -block.blockId, Edge::BLACK);
			}
		}

		//chromosome ends
		Edge* headEdge = this->addEdge(INFINUM, perm.blocks.front().signedId(),
									   perm.seqId);
		headEdge->rightPos = perm.blocks.front().start;
		Edge* tailEdge = this->addEdge(-perm.blocks.back().signedId(), INFINUM,
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

			maxBlockOverlap = std::max(maxBlockOverlap, leftPos - rightPos);
			if (rightPos < leftPos) leftPos = rightPos;

			/*if (rightPos < leftPos)
				std::cerr << "WARNING: overlapping blocks\n"
						  << blockPair.first.start << " " << blockPair.first.end
						  << " " << blockPair.second.start << " "
						  << blockPair.second.end << " | " << perm.seqId << "\n";*/

			curEdge = this->addEdge(-blockPair.first.signedId(),
									blockPair.second.signedId(), perm.seqId);
			curEdge->leftPos = leftPos;
			curEdge->rightPos = rightPos;

			prevEdge->nextEdge = curEdge;
			curEdge->prevEdge = prevEdge;
			prevEdge = curEdge;
		}

		curEdge->nextEdge = tailEdge;
		tailEdge->prevEdge = curEdge;
	}

	if (maxBlockOverlap > 0)
	{
		std::cerr << "\n\tWARNING: some alignment blocks were overlapping "
			"by at most " << maxBlockOverlap << "\n" <<
			"\tIt is unexpected and might result into nonsense synteny blocks.\n" <<
			"\tPlease note that this tool currently support alignments\n" <<
			"\tproduced by either Cactus or SibeliaZ.\n\n";

	}
	DEBUG_PRINT("Constructed graph with " << _nodes.size() << " nodes");
}

/*namespace
{
	std::unordered_map<Edge*, int> getConjunctionEdges(BreakpointGraph& bg)
	{
		DEBUG_PRINT("Getting conjunction edges");
		std::unordered_map<Edge*, int> edgeToGroup;
		std::unordered_map<Edge*, SetNode<int>*> setNodes;
		int nextId = 1;
		auto getSetNode = [&setNodes, &nextId] (Edge* e) 
		{
			if (!setNodes.count(e))
			{
				setNodes[e] = new SetNode<int>(nextId++);
			}
			return setNodes[e];
		};

		for (int node : bg.iterNodes())
		{
			if (!bg.isBifurcation(node)) continue;

			NodeVec neighbors = bg.getNeighbors(node);
			if (neighbors.size() != 3 || !contains(neighbors, bg.INFINUM))
			{
				continue;
			}
			neighbors.erase(std::remove(neighbors.begin(), neighbors.end(), 
										int(bg.INFINUM)), neighbors.end());
			
			Edge* leftEdge = bg.getAdjacentBlackEdge(neighbors[0]);
			Edge* rightEdge = bg.getAdjacentBlackEdge(neighbors[1]);
			unionSet(getSetNode(leftEdge), getSetNode(rightEdge));
		}

		for (auto &nodePair : setNodes)
		{
			edgeToGroup[nodePair.first] = findSet(nodePair.second)->data;
		}
		for (auto &nodePair : setNodes) delete nodePair.second;

		DEBUG_PRINT("Getting conjunction edges - done");
		return edgeToGroup;
	}
}*/

void BreakpointGraph::getPermutations(PermVec& permutations, 
									  BlockGroups& blockGroups)
{
	DEBUG_PRINT("Reading permutations from graph");

	int nextEdgeId = 1;
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
		Edge* curEdge = startEdge;
		do
		{
			Edge* prevEdge = curEdge;
			curEdge = curEdge->nextEdge;
			EdgeVec blackEdges = this->getBlackEdges(prevEdge->rightNode, 
													 curEdge->leftNode);
			assert(blackEdges.size() == 1);

			int blockId = getEdgeId(blackEdges[0]);
			int start = prevEdge->rightPos;
			int end = curEdge->leftPos;
			assert(end >= start);
			int sign = (blackEdges[0]->rightNode == curEdge->leftNode) ? 1 : -1;

			permutations.back().blocks.push_back(Block(blockId, sign, 
													   start, end));
			permutations.back().seqId = startEdge->seqId;
			permutations.back().seqName = _seqNames[startEdge->seqId];
			permutations.back().nucLength = _seqLength[startEdge->seqId];
		}
		while (!curEdge->hasNode(INFINUM));
	}
	//std::unordered_map<Edge*, int> edgeToGroup = getConjunctionEdges(*this);
	//for (auto &edgePair : edgeToGroup)
	//	blockGroups[getEdgeId(edgePair.first)] = edgePair.second;
		
	DEBUG_PRINT("Reading permutations from graph - finished");
}
