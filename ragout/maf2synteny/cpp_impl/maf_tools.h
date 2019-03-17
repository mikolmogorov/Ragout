//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include "breakpoint_graph.h"

PermVec mafToPermutations(const std::string& mafFile, int minBlockLen);

PermVec parseGff(const std::string& filename, int minBlockLen);
