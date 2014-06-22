//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>
#include <list>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <cassert>
#include <iostream>

#include "permutation.h"
#include "utility.h"

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
	size_t leftPos;
	size_t rightPos;
	Edge* prevEdge;
	Edge* nextEdge;
};

typedef std::vector<Edge*> 	EdgeVec;
typedef std::vector<int> 	NodeVec;

struct Node
{
	EdgeVec edges;
	NodeVec neighbors;
};

class BreakpointGraph
{
public:
	BreakpointGraph(const std::vector<Permutation>& permutations);
	~BreakpointGraph();

	Edge* addEdge(int leftNode, int rightNode, int seqId)
	{
		assert((leftNode != rightNode) || seqId != Edge::BLACK);
		Edge* edge = new Edge(leftNode, rightNode, seqId);

		_nodes[leftNode].edges.push_back(edge);
		if (!contains(_nodes[leftNode].neighbors, rightNode))
		{
			_nodes[leftNode].neighbors.push_back(rightNode);
		}

		if (leftNode != rightNode)
		{
			_nodes[rightNode].edges.push_back(edge);
			if (!contains(_nodes[rightNode].neighbors, leftNode))
			{
				_nodes[rightNode].neighbors.push_back(leftNode);
			}
		}
		
		return edge;
	}

	EdgeVec getEdges(int nodeOne, int nodeTwo)
	{
		assert(nodeOne && nodeTwo);
		EdgeVec edgesOut;
		for (Edge* e : _nodes[nodeOne].edges)
		{
			if (e->hasNode(nodeTwo)) edgesOut.push_back(e);
		}
		return edgesOut;
	}

	EdgeVec getBlackEdges(int nodeOne, int nodeTwo)
	{
		assert(nodeOne && nodeTwo);
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
		assert(nodeOne && nodeTwo);
		EdgeVec edgesOut;
		for (Edge* e : _nodes[nodeOne].edges)
		{
			if (e->hasNode(nodeTwo) && e->seqId != Edge::BLACK)
				edgesOut.push_back(e);
		}
		return edgesOut;
	}

	void removeEdges(int idOne, int idTwo)
	{
		//assert(idOne != idTwo);

		Node& nodeOne = _nodes.at(idOne);
		Node& nodeTwo = _nodes.at(idTwo);

		assert(contains(nodeOne.neighbors, idTwo));
		assert(contains(nodeTwo.neighbors, idOne));
		vecRemove(nodeOne.neighbors, idTwo);
		vecRemove(nodeTwo.neighbors, idOne);

		auto e = nodeOne.edges.begin();
		while (e != nodeOne.edges.end()) 
		{
			if ((*e)->hasNode(idTwo))
			{
				delete *e;
				if (idOne != idTwo) vecRemove(nodeTwo.edges, *e);
				e = nodeOne.edges.erase(e);
			}
			else
			{
				++e;
			}
		}
	}

	void removeNode(int node)
	{
		if (!_nodes.count(node))
			return;
		for (int neighbor : this->getNeighbors(node))
			this->removeEdges(neighbor, node);
		_nodes.erase(node);
	}

	Edge* getAdjacentBlackEdge(int node)
	{
		for (Edge* e : _nodes[node].edges)
		{
			if (e->seqId == Edge::BLACK)
				return e;
		}
		return nullptr;
	}

	const NodeVec& getNeighbors(int node)
	{
		assert(node);
		return _nodes[node].neighbors;
	}

	bool isBifurcation(int node)
	{
		NodeVec neighbors = this->getNeighbors(node);
		if (neighbors.size() > 2)
			return true;

		//all edges should be either balck or colred
		for (auto neighbor : neighbors)
		{
			EdgeVec edges = this->getEdges(node, neighbor);

			bool hasBlack = false;
			for (auto e : edges)
				if (e->seqId == Edge::BLACK)
					hasBlack = true;
			if (hasBlack && edges.size() > 1)
				return true;
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

	void getPermutations(PermVec& permutations, BlockGroups& blockGroups);

	static const int INFINUM = std::numeric_limits<int>::max();

private:
	std::unordered_map<int, Node> 	_nodes;
	std::vector<Edge*> 				_origins;
	std::unordered_map<int, int> 	_seqLength;
	std::unordered_map<int, std::string> _seqNames;

};
