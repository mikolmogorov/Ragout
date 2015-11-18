//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#include "breakpoint_graph.h"
#include "compress_algorithms.h"
#include "utility.h"

#include <deque>
#include <unordered_set>
#include <iostream>

bool extendPath(BreakpointGraph& graph, int prevNode, int curNode, 
				int maxGap, std::deque<int>& outPath)
{
	outPath.push_back(prevNode);
	outPath.push_back(curNode);
	while(true)
	{
		//check for bifurcation
		if (graph.isBifurcation(curNode) || graph.INFINUM == prevNode ||
			graph.INFINUM == curNode)
			return true;

		//check distance
		for (Edge* e : graph.getColoredEdges(prevNode, curNode))
			if (e->getLen() > maxGap)
				return false;

		//everything is ok, continue path
		NodeVec neighbors = graph.getNeighbors(curNode);
		assert(neighbors.size() == 2);
		int otherNode = (neighbors[0] != prevNode) ? neighbors[0] : neighbors[1];
		
		prevNode = curNode;
		curNode = otherNode;
		outPath.push_back(curNode);
	}
}

//The code below was written in military registration office
//during my waiting for departure to Scientific Detachment

bool compressPath(BreakpointGraph& graph, std::deque<int>& path,
				  std::unordered_set<int>& nodesToDel)
{
	if (path.size() <= 2) return false;

	//ensure we start and end with black edges
	if (graph.getBlackEdges(path[0], path[1]).empty())
		path.pop_front();
	if (graph.getBlackEdges(path[path.size()-2], path[path.size()-1]).empty())
		path.pop_back();
	if (path.size() == 2) return false;
	assert (path.size() >= 4);
	
	//updating graph
	graph.removeEdges(path[0], path[1]);
	graph.removeEdges(path[path.size()-2], path[path.size()-1]);
	std::copy(path.begin() + 1, path.end() - 1, 
			  std::inserter(nodesToDel, nodesToDel.begin()));
	graph.addEdge(path.front(), path.back(), Edge::BLACK);

	//updating links
	//for each colored edge in path
	for (Edge* adj : graph.getColoredEdges(path[1], path[2]))
	{
		Edge* headAdj = adj->nextEdge;
		while(!headAdj->hasNode(path.front()) && !headAdj->hasNode(path.back()))
			headAdj = headAdj->nextEdge;

		Edge* tailAdj = adj->prevEdge;
		while(!tailAdj->hasNode(path.front()) && !tailAdj->hasNode(path.back()))
			tailAdj = tailAdj->prevEdge;

		assert(headAdj->seqId == tailAdj->seqId);
		headAdj->prevEdge = tailAdj;
		tailAdj->nextEdge = headAdj;
	}

	return true;
}

int compressGraph(BreakpointGraph& graph, int maxGap)
{
	DEBUG_PRINT("Started compression");
	int numCompressed = 0;
	std::unordered_set<int> nodesToDel;

	int failsLength = 0;
	int failsStructure = 0;

	for (int node : graph.iterNodes())
	{
		if (nodesToDel.count(node) || !graph.isBifurcation(node))
			continue;

		for (int neighbor : graph.getNeighbors(node))
		{
			assert(graph.getNeighbors(node).size() >= 2);
			std::deque<int> path;
			if (extendPath(graph, node, neighbor, maxGap, path))
			{
				failsStructure += 1;
			}
			else
			{
				failsLength += 1;
			}
			if (compressPath(graph, path, nodesToDel))
			{
				++numCompressed;
			}
		}
	}

	//cleaning up
	for (int node : nodesToDel)
	{
		graph.removeNode(node);
	}
	
	DEBUG_PRINT("Stuctural fails: " << failsStructure);
	DEBUG_PRINT("Length fails: " << failsLength);
	DEBUG_PRINT("Finished compression");
	return numCompressed;
}

typedef std::vector<std::deque<int>> BranchSet;

bool collapseBulge(BreakpointGraph& graph, const BranchSet& branches,
				   std::unordered_set<int>& nodesToDel)
{
	//checking bulge structure
	for (auto &branch : branches)
	{
		if (branch.size() != 2 && branch.size() != 4)
			return false;
	}

	for (auto &branch : branches)
	{
		if (branch.size() == 2) //a branch from one colored edge, nothing to do
			continue;

		//replace three edges with one (colored) -- "hiding" varinance
		for (Edge* adj : graph.getColoredEdges(branch[0], branch[1]))
		{
			Edge* nextAdj = adj->nextEdge;
			Edge* newAdj = nullptr;
			Edge* prevAdj = nullptr;
			//checking direction
			if (contains(branch, nextAdj->leftNode))
			{
				nextAdj = nextAdj->nextEdge;
				prevAdj = adj->prevEdge;

				newAdj = graph.addEdge(branch.front(), branch.back(), 
									   adj->seqId);
				newAdj->leftPos = adj->leftPos;
				newAdj->rightPos = adj->nextEdge->rightPos;
			}
			else
			{
				prevAdj = adj->prevEdge->prevEdge;

				newAdj = graph.addEdge(branch.back(), branch.front(),
									   adj->seqId);
				newAdj->leftPos = adj->prevEdge->leftPos;
				newAdj->rightPos = adj->rightPos;
			}

			newAdj->nextEdge = nextAdj;
			newAdj->prevEdge = prevAdj;
			nextAdj->prevEdge = newAdj;
			prevAdj->nextEdge = newAdj;
		}

		//cleaning up
		assert(!graph.getEdges(branch[0], branch[1]).empty());
		graph.removeEdges(branch[0], branch[1]);
		graph.removeEdges(branch[2], branch[3]);

		std::copy(branch.begin() + 1, branch.end() - 1,
				  std::inserter(nodesToDel, nodesToDel.begin()));
	}

	return true;
}

bool findBulge(BreakpointGraph& graph, int node, int maxGap, BranchSet& branches)
{
	std::unordered_map<int, BranchSet> byEnd;
	for (int neighbor : graph.getNeighbors(node))
	{
		std::deque<int> path;
		extendPath(graph, node, neighbor, maxGap, path);
		byEnd[path.back()].push_back(std::move(path));
	}

	//checking bulge structure: -<=>-
	if (byEnd.size() != 2) return false;

	EdgeVec leftFlank;
	EdgeVec rightFlank;
	for (auto endPair : byEnd)
	{
		if (endPair.second.size() == 1)
			leftFlank = graph.getBlackEdges(endPair.second[0][0], 
											endPair.second[0][1]);
		else
			branches = endPair.second;
	}
	if (branches.empty() || leftFlank.empty()) return false;

	int pathEnd = branches.front().back();
	if (pathEnd == node) return false;

	std::vector<int> branchRepr;
	for (auto branch : branches)
		branchRepr.push_back(branch[branch.size()-2]);
	std::vector<int> otherNeighbors;
	for (int neighbor : graph.getNeighbors(pathEnd))
	{
		if (!contains(branchRepr, neighbor))
			otherNeighbors.push_back(neighbor);
	}
	if (otherNeighbors.size() > 1 || otherNeighbors.empty())
		return false;

	rightFlank = graph.getBlackEdges(pathEnd, otherNeighbors.front());
	if (rightFlank.empty())
		return false;

	return true;
}

int removeBulges(BreakpointGraph& graph, int maxGap)
{
	DEBUG_PRINT("Removing bulges");
	int numCollapsed = 0;
	std::unordered_set<int> nodesToDel;
	for (int node : graph.iterNodes())
	{
		if (nodesToDel.count(node) || !graph.isBifurcation(node))
			continue;

		BranchSet branches;
		if (findBulge(graph, node, maxGap, branches))
			if (collapseBulge(graph, branches, nodesToDel))
				++numCollapsed;
	}
	for (int node : nodesToDel)
	{
		graph.removeNode(node);
	}

	DEBUG_PRINT("Removing bulges - done");
	return numCollapsed;
}
