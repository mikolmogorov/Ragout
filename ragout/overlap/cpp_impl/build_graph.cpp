//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

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
#include "overlap.h"

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

bool getContigs(const std::string& filename, std::vector<FastaRecord>& contigs)
{
	std::cerr << "\tReading FASTA\n";
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
	std::cerr << "\tKmer size is set to " << mostFrequent << std::endl;

	std::copy_if(overlapsIn.begin(), overlapsIn.end(),
				 std::back_inserter(overlapsOut),
				 [&] (const Overlap& o) {return o.size == mostFrequent;});
	return overlapsOut;
}

void drawHistogram(std::vector<Overlap>& overlaps)
{
	std::map<int, int> hist;
	for (Overlap& ovlp : overlaps)
	{
		++hist[ovlp.size];
	}
	for (auto histPair : hist)
	{
		std::cerr << histPair.first << " " << histPair.second << std::endl;
	}
}

}	//end anonymous namespace

bool makeOverlapGraph(const std::string& fileIn, const std::string& fileOut, 
		  			  int minOverlap, int maxOverlap, bool filterKmer,
					  bool drawHist)
{
	std::vector<FastaRecord> contigs;
	std::ofstream streamOut(fileOut);
	if (!streamOut)
	{
		std::cerr << "Cannot open: " << fileOut << std::endl;
		return false;
	}
	if (!getContigs(fileIn, contigs)) return false;

	std::vector<Overlap> overlaps = getOverlaps(contigs, minOverlap, 
												maxOverlap);
	if (drawHist)
	{
		drawHistogram(overlaps);
	}
	if (filterKmer)
	{
		overlaps = filterByKmer(overlaps);
	}

	buildAssemblyGraph(contigs, overlaps, streamOut);
	return true;
}
