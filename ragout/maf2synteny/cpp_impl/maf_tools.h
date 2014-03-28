#pragma once

#include <string>
#include "breakpoint_graph.h"

std::vector<Permutation> mafToPermutations(const std::string& mafFile, 
										   int minBlockLen);
