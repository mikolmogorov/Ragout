#pragma once

#include "breakpoint_graph.h"

int compressGraph(BreakpointGraph& graph, int maxGap);
int removeBulges(BreakpointGraph& graph, int maxGap);
