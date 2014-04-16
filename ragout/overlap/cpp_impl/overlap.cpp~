#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <utility>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <map>

#include "fasta.h"

namespace
{

struct Edge
{
	Edge(int begin, int end, const std::string& contigId):
		begin(begin), end(end), contigId(contigId) {}
	int begin;
	int end;
	std::string contigId;
};

void makePrefixFunction(const std::string& str, std::vector<int>& table)
{
	table.assign(str.length(), 0);
	for (size_t i = 1; i < str.length(); ++i)
	{
		int k = table[i - 1];
		while (true)
		{
			if (str[i] == str[k])
			{
				table[i] = k + 1;
				break;
			}
			if (k == 0)
			{
				table[i] = 0;
				break;
			}
			k = table[k - 1];
		}
	}
}


int kmpOverlap(const std::string& strHead, const std::string& strTail, 
			   const std::vector<int>& prefixFun)
{
	size_t i = 0;
	size_t j = 0;
	size_t maxOvlp = 0;
	while (i + j < strHead.length())
	{
		if (strHead[i + j] == strTail[j])
		{
			if (i + j == strHead.length() - 1) 
			{
				maxOvlp = std::max(maxOvlp, j + 1);
			}
			
			if (j == strTail.length() - 1)
			{
				i = i + j - prefixFun[j - 1];
				j = prefixFun[j - 1];
			}
			else
			{
				++j;	
			}
		}
		else if (j == 0)
		{
			++i;
		}
		else
		{
			i = i + j - prefixFun[j - 1];
			j = prefixFun[j - 1];
		}
	}
	return maxOvlp;
}

int getOverlap(FastaRecord* headContig, FastaRecord* tailContig, size_t maxOvlp)
{
	static std::unordered_map<FastaRecord*, std::vector<int>> prefixFunCache;

	maxOvlp = std::min(maxOvlp, headContig->sequence_.length());
	maxOvlp = std::min(maxOvlp, tailContig->sequence_.length());
	std::string strHead(headContig->sequence_, 
						headContig->sequence_.length() - maxOvlp);
	std::string strTail(tailContig->sequence_, 0, maxOvlp);

	auto &pfun = prefixFunCache[tailContig];
	if (pfun.size() != maxOvlp)	//zero, for instance
	{
		makePrefixFunction(strTail, pfun);
	}
	return kmpOverlap(strHead, strTail, pfun);
}

int newNodeId()
{
	static int nodeId = 0;
	return nodeId++;
}

bool getContigs(const std::string& filename, std::vector<FastaRecord>& contigs)
{
	FastaReader reader(filename);
	if (!reader.IsOk()) 
	{
		std::cerr << "Error openning " << filename << std::endl;
		return false;
	}
	try
	{
		reader.GetSequencesWithComplements(contigs);
	}
	catch (std::runtime_error& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return false;
	}
	return true;
}

bool buildGraph(const std::string& filename, int minOverlap, 
				int maxOverlap, std::list<Edge>& edges)
{
	//std::map<int, int> ovlpHist;

	std::vector<FastaRecord> contigs;
	if (!getContigs(filename, contigs)) return false;
	
	std::unordered_map<FastaRecord*, int> contigNodes;
	std::unordered_set<FastaRecord*> visited;
	std::vector<std::pair<FastaRecord*, int>> dfsStack;

	int prevProgress = 0;
	std::cerr << "Progress: ";

	for (auto &contig : contigs)
	{
		if (visited.count(&contig)) continue;

		dfsStack.push_back(std::make_pair(&contig, newNodeId()));
		visited.insert(&contig);
		while(!dfsStack.empty())
		{
			FastaRecord* curContig = dfsStack.back().first;
			int leftNodeId = dfsStack.back().second;
			dfsStack.pop_back();
			contigNodes[curContig] = leftNodeId;

			//finding overlaps
			std::list<FastaRecord*> overlaps;
			for (auto &otherContig : contigs)
			{
				if (curContig == &otherContig) continue;

				int overlap = getOverlap(curContig, &otherContig, maxOverlap);
				if (overlap >= minOverlap)
				{
					overlaps.push_back(&otherContig);
					//ovlpHist[overlap] += 1;
				}
			}

			//processing them
			int rightNodeId = -1;
			if (!overlaps.empty())
			{
				FastaRecord* sampleContig = overlaps.front();
				auto foundNode = contigNodes.find(sampleContig);
				if (foundNode != contigNodes.end())
				{
					rightNodeId = foundNode->second;
				}
				else
				{
					rightNodeId = newNodeId();
					for (auto ovlp : overlaps) contigNodes[ovlp] = rightNodeId;
				}

				for (auto ovlp : overlaps)
				{
					if (visited.find(ovlp) == visited.end())
					{
						dfsStack.push_back(std::make_pair(ovlp, rightNodeId));
						visited.insert(ovlp);
					}
				}
			}
			else
			{
				rightNodeId = newNodeId();
			}

			edges.push_back(Edge(leftNodeId, rightNodeId, curContig->description_));

			int progress = edges.size() * 100 / contigs.size();
			if (progress / 10 > prevProgress / 10)
			{
				std::cerr << progress << " ";
				prevProgress = progress;
			}
		}
	}
	std::cerr << std::endl;
	//for (auto histPair : ovlpHist)
	//{
	//	std::cerr << histPair.first << " " << histPair.second << std::endl;
	//}
	return true;
}

}	//end anonymous namespace

bool makeOverlapGraph(const std::string& fileIn, const std::string& fileOut, 
		  			  int minOverlap, int maxOverlap)
{
	std::list<Edge> edges;
	if (!buildGraph(fileIn, minOverlap, maxOverlap, edges)) return false;
	
	std::ofstream streamOut(fileOut);
	if (!streamOut)
	{
		std::cerr << "Cannot open: " << fileOut << std::endl;
		return false;
	}

	streamOut << "digraph {\n";
	for (auto edge : edges)
	{
		streamOut << edge.begin << " -> " << edge.end 
				  << " [label=\"" << edge.contigId << "\"];\n";
	}
	streamOut << "}\n";

	return true;
}


