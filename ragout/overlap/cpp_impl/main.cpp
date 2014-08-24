//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <vector>
#include <string>
#include <string.h>

#include "build_graph.h"

void printUsage()
{
	std::cerr << "Usage: overlap fasta_in dot_out min_k "
			  << "max_k [--detect-kmer] [--hist]\n"
	 		  << "Constructs overlap graph from input contigs"
			  << "and outputs it in dot format\n";
}

int main(int argc, char** argv)
{
	std::vector<std::string> posArgs;
	bool detectKmer = false;
	bool drawHist = false;
	for (char** arg = argv; arg != argv + argc; ++arg)
	{
		if (strcmp(*arg, "--detect-kmer") == 0)
		{
			detectKmer = true;
		}
		else if (strcmp(*arg, "--hist") == 0)
		{
			drawHist = true;
		}
		else if (strcmp(*arg, "--help") == 0)
		{
			printUsage();
			return 0;
		}
		else
		{
			posArgs.push_back(*arg);
		}
	}

	if (argc < 5)
	{
		printUsage();
		return 1;
	}

	return !makeOverlapGraph(posArgs[1], posArgs[2], atoi(posArgs[3].c_str()),
							 atoi(posArgs[4].c_str()), detectKmer, drawHist);
}
