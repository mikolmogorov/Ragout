//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <vector>
#include <string>

#include "build_graph.h"

void printUsage()
{
	std::cerr << "Usage: overlap fasta_in dot_out min_k "
			  << "max_k [--detect-kmer]\n"
	 		  << "Constructs overlap graph from input contigs"
			  << "and outputs it in dot format\n";
}

int main(int argc, char** argv)
{
	std::vector<std::string> args(argv, argv + argc);
	bool detectKmer = false;
	for (auto itArg = args.begin(); itArg != args.end(); ++itArg)
	{
		if (*itArg == "--detect-kmer")
		{
			detectKmer = true;
			args.erase(itArg);
			break;
		}
		if (*itArg == "--help")
		{
			printUsage();
			return 0;
		}
	}

	if (argc < 5)
	{
		printUsage();
		return 1;
	}

	return !makeOverlapGraph(args[1], args[2], atoi(args[3].c_str()),
							 atoi(args[4].c_str()), detectKmer);
}
