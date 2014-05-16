#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <stdexcept>

#include "fasta.h"
#include "disjoint_set.h"

namespace
{

struct Edge
{
	Edge(const std::string& begin, const std::string& end, 
		 const std::string& label): begin(begin), end(end), label(label) {}
	std::string begin;
	std::string end;
	std::string label;
};

struct Overlap
{
	Overlap(FastaRecord* prevContig, FastaRecord* nextContig, int size):
		prevContig(prevContig), nextContig(nextContig), size(size) {}

	FastaRecord* prevContig;
	FastaRecord* nextContig;
	int size;
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


std::vector<Overlap> overlapAll(std::vector<FastaRecord>& contigs, 
								int minOverlap, int maxOverlap)
{
	std::vector<Overlap> overlaps;

	int prevProgress = 0;
	int contigCounter = 0;
	std::map<int, int> overlapHist;
	std::cerr << "Progress: ";

	for (auto &headContig : contigs)
	{
		for (auto &tailContig : contigs)
		{
			if (&headContig == &tailContig) continue;
			
			int overlapLen = getOverlap(&headContig, &tailContig, maxOverlap);
			if (overlapLen >= minOverlap)
			{
				overlaps.push_back(Overlap(&headContig, &tailContig, 
										   overlapLen));
				overlapHist[overlapLen] += 1;
			}
		}
		
		++contigCounter;
		int progress = contigCounter * 100 / contigs.size();
		if (progress / 10 > prevProgress / 10)
		{
			std::cerr << progress << " ";
			prevProgress = progress;
		}
	}
	std::cerr << std::endl;
	for (auto ovlpPair : overlapHist)
	{
		std::cerr << ovlpPair.first << " " << ovlpPair.second << std::endl;
	}

	return overlaps;
}


void buildAssemblyGraph(std::vector<FastaRecord>& contigs,
						std::vector<Overlap>& overlaps, std::ostream& streamOut)
{
	//std::unordered_set<FastaRecord*> usedContigs;

	streamOut << "digraph {\n";
	for (Overlap& ovlp : overlaps)
	{
		//usedContigs.insert(ovlp.prevContig);
		//usedContigs.insert(ovlp.nextContig);

		streamOut << "\"" << ovlp.prevContig->description_ << "\" -> \""
				  << ovlp.nextContig->description_ 
				  << "\" [label=\"" << ovlp.size << "\"];\n";
	}

	//for (FastaRecord& rec : contigs)
	//{
	//	if (!usedContigs.count(&rec))
	//	{
	//		streamOut << "\"" << rec.description_ << "\"\n";
	//	}
	//}

	streamOut << "}\n";
}

std::vector<Overlap> filterByKmer(std::vector<Overlap>& overlapsIn)
{
	std::vector<Overlap> overlapsOut;
	std::unordered_map<int, int> overlapHist;
	
	int mostFrequent = -1;
	int frequency = 0;
	for (auto &ovlp : overlapsIn) ++overlapHist[ovlp.size];
	for (auto &histPair : overlapHist)
	{
		if (histPair.second > frequency)
		{
			mostFrequent = histPair.first;
			frequency = histPair.second;
		}
	}
	std::cerr << "Kmer size is set to " << mostFrequent << std::endl;

	std::copy_if(overlapsIn.begin(), overlapsIn.end(),
				 std::back_inserter(overlapsOut),
				 [&] (const Overlap& o) {return o.size == mostFrequent;});
	return overlapsOut;
}

}	//end anonymous namespace

bool makeOverlapGraph(const std::string& fileIn, const std::string& fileOut, 
		  			  int minOverlap, int maxOverlap, bool filterKmer)
{
	std::vector<FastaRecord> contigs;
	std::ofstream streamOut(fileOut);
	if (!streamOut)
	{
		std::cerr << "Cannot open: " << fileOut << std::endl;
		return false;
	}
	if (!getContigs(fileIn, contigs)) return false;

	std::vector<Overlap> overlaps = overlapAll(contigs, minOverlap, maxOverlap);
	if (filterKmer)
	{
		overlaps = filterByKmer(overlaps);
	}

	buildAssemblyGraph(contigs, overlaps, streamOut);
	return true;
}
