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
	//Node(int id): nodeId(id) {}
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

	void getEdges(int nodeOne, int nodeTwo, std::vector<Edge*>& edgesOut)
	{
		for (Edge* e : _nodes[nodeOne].edges)
		{
			if (e->hasNode(nodeTwo)) edgesOut.push_back(e);
		}
	}

	void getBlackEdges(int nodeOne, int nodeTwo, std::vector<Edge*>& edgesOut)
	{
		for (Edge* e : _nodes[nodeOne].edges)
		{
			if (e->hasNode(nodeTwo) && e->seqId == Edge::BLACK)
				edgesOut.push_back(e);
		}
	}

	void getColoredEdges(int nodeOne, int nodeTwo, std::vector<Edge*>& edgesOut)
	{
		for (Edge* e : _nodes[nodeOne].edges)
		{
			if (e->hasNode(nodeTwo) && e->seqId != Edge::BLACK)
				edgesOut.push_back(e);
		}
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

	void getNeighbors(int node, std::vector<int>& outNodes)
	{
		std::unordered_set<int> neighbors;
		for (auto e : _nodes[node].edges)
		{
			neighbors.insert(e->leftNode != node ? e->leftNode : e->rightNode);
		}
		std::copy(neighbors.begin(), neighbors.end(), 
				  std::back_inserter(outNodes));
	}

	bool isBifurcation(int node)
	{
		std::vector<int> neighbors;
		this->getNeighbors(node, neighbors);
		if (neighbors.size() > 2) return true;

		//all edges should be either balck or colred
		for (auto neighbor : neighbors)
		{
			std::vector<Edge*> edges;
			this->getEdges(node, neighbor, edges);

			bool hasBlack = false;
			for (auto e : edges) if (e->seqId == Edge::BLACK) hasBlack = true;
			if (hasBlack && edges.size() > 1) return true;
		}
		return false;
	}

	void getFragmentedBlocks(std::vector<std::vector<int>>& groups);
	void getPermutations(std::vector<Permutation>& permutstions);

	static const int INFINUM = std::numeric_limits<int>::max();

private:
	std::unordered_map<int, Node> 	_nodes;
	std::vector<Edge*> 				_origins;

};
