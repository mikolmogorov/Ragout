//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "breakpoint_graph.h"

int compressGraph(BreakpointGraph& graph, int maxGap);
int removeBulges(BreakpointGraph& graph, int maxGap);
