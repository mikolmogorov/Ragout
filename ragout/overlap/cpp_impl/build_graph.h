//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#ifndef BUILD_GRAPH_H
#define BUILD_GRAPH_H

#include <string>

bool makeOverlapGraph(const std::string& fileIn, const std::string& fileOut, 
		  			  int minOverlap, int maxOverlap, bool filterKmer,
					  bool drawHist);

#endif
