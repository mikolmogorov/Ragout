#pragma once

#include <string>
#include "breakpoint_graph.h"

PermVec mafToPermutations(const std::string& mafFile, int minBlockLen);
