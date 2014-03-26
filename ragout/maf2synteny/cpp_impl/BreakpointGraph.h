#pragma once

#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <cassert>

struct Edge
{
	Edge(int left, int right, int seqId): 
		leftNode(left), rightNode(right), seqId(seqId),
		prevEdge(nullptr), nextEdge(nullptr) {}
	bool hasNode(int node) 	{return node == leftNode || node == rightNode;}
	int  getLen() 			{return rightPos - leftPos;}

	static const int BLACK = -1;

	int leftNode;
	int rightNode;
	int seqId;
	int leftPos;
	int rightPos;
	Edge* prevEdge;
	Edge* nextEdge;
};

struct Node
{
	std::list<Edge*> edges;
};

struct Block
{
	Block(int blockId, int sign, int start, int end):
		blockId(blockId), sign(sign), start(start), end(end) {}
	int getLen() const {assert(end >= start); return end - start;}

	int blockId;
	int sign;
	int start;
	int end;
};

struct Permutation
{
	int seqId;
	std::vector<Block> blocks;
};

typedef std::vector<Edge*> EdgeVec;
typedef std::vector<int> NodeVec;

class BreakpointGraph
{
public:
	BreakpointGraph(const std::vector<Permutation>& permutations);

	Edge* addEdge(int leftNode, int rightNode, int seqId)
	{
		Edge* edge = new Edge(leftNode, rightNode, seqId);
		_nodes[leftNode].edges.push_back(edge);
		_nodes[rightNode].edges.push_back(edge);
		return edge;
	}

	EdgeVec getEdges(int nodeOne, int nodeTwo)
	{
		EdgeVec edgesOut;
		for (Edge* e : _nodes[nodeOne].edges)
		{
			if (e->hasNode(nodeTwo)) edgesOut.push_back(e);
		}
		return edgesOut;
	}

	EdgeVec getBlackEdges(int nodeOne, int nodeTwo)
	{
		EdgeVec edgesOut;
		for (Edge* e : _nodes[nodeOne].edges)
		{
			if (e->hasNode(nodeTwo) && e->seqId == Edge::BLACK)
				edgesOut.push_back(e);
		}
		return edgesOut;
	}

	EdgeVec getColoredEdges(int nodeOne, int nodeTwo)
	{
		EdgeVec edgesOut;
		for (Edge* e : _nodes[nodeOne].edges)
		{
			if (e->hasNode(nodeTwo) && e->seqId != Edge::BLACK)
				edgesOut.push_back(e);
		}
		return edgesOut;
	}

	void removeEdges(int nodeOne, int nodeTwo)
	{
		for (auto e = _nodes[nodeOne].edges.begin();
			 	e != _nodes[nodeOne].edges.end(); ++e)
		{
			if ((*e)->hasNode(nodeTwo))
			{
				e = _nodes[nodeOne].edges.erase(e);
				_nodes[nodeTwo].edges.remove(*e);
			}
		}
	}

	NodeVec getNeighbors(int node)
	{
		NodeVec outNodes;
		std::unordered_set<int> neighbors;
		for (auto e : _nodes[node].edges)
		{
			neighbors.insert(e->leftNode != node ? e->leftNode : e->rightNode);
		}
		std::copy(neighbors.begin(), neighbors.end(), 
				  std::back_inserter(outNodes));
		return outNodes;
	}

	bool isBifurcation(int node)
	{
		NodeVec neighbors = this->getNeighbors(node);
		if (neighbors.size() > 2) return true;

		//all edges should be either balck or colred
		for (auto neighbor : neighbors)
		{
			EdgeVec edges = this->getEdges(node, neighbor);

			bool hasBlack = false;
			for (auto e : edges) if (e->seqId == Edge::BLACK) hasBlack = true;
			if (hasBlack && edges.size() > 1) return true;
		}
		return false;
	}

	NodeVec iterNodes()
	{
		NodeVec nodes;
		for (auto itNode : _nodes)
			nodes.push_back(itNode.first);
		return nodes;
	}

	void getFragmentedBlocks(std::vector<std::vector<int>>& groups);
	void getPermutations(std::vector<Permutation>& permutstions);

	static const int INFINUM = std::numeric_limits<int>::max();

private:
	std::unordered_map<int, Node> 	_nodes;
	std::vector<Edge*> 				_origins;

};
