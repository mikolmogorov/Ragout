#include "BreakpointGraph.h"

void extendPath(BreakpointGraph& graph, int prevNode, int curNode, 
				int maxGap, std::vector<int>& outPath)
{
	outPath.push_back(prevNode);
	outPath.push_back(curNode);
	while(true)
	{
		//chack for bifurcation
		if (graph.isBifurcation(curNode) || graph.INFINUM == prevNode ||
			graph.INFINUM == curNode)
			break;

		//check distance
		std::vector<Edge*> edges;
		graph.getColoredEdges(prevNode, curNode, edges);
		for (auto e : edges)
			if (e.getLen() > maxGap)
				break;

		//everything is ok, continue path
		std::vector<int> neighbors;
		graph.getNeighbors(curNode);
		int otherNode = (neighbors[0] == prevNode) ? neighbors[0] : neighbors[1];
		
		outPath.push_back(curNode);
		prevNode = curNode;
		curNode = otherNode;
	}
}
