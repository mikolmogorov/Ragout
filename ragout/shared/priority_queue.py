#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
A priority queue implementation.
Based on http://stackoverflow.com/questions/407734/a-generic-priority-queue-for-python
"""
import heapq

class PriorityQueue(object):
    """
    Priority queue based on heap, capable of inserting a new node with
    desired priority, updating the priority of an existing node and deleting
    an abitrary node while keeping invariant
    """

    def __init__(self):
        """
        if 'heap' is not empty, make sure it's heapified
        """
        self.heap = []
        self.entry_finder = {}
        self.REMOVED = "<remove_marker>"
        self.length = 0

    def get_length(self):
        return self.length

    def insert(self, node, priority=0):
        """
        'entry_finder' bookkeeps all valid entries, which are bonded in
        'heap'. Changing an entry in either leads to changes in both.
        """
        if node in self.entry_finder:
            self.delete(node)
        entry = [priority, node]
        self.entry_finder[node] = entry
        heapq.heappush(self.heap, entry)
        self.length += 1

    def delete(self, node):
        """
        Instead of breaking invariant by direct removal of an entry, mark
        the entry as "REMOVED" in 'heap' and remove it from 'entry_finder'.
        Logic in 'pop()' properly takes care of the deleted nodes.
        """
        entry = self.entry_finder.pop(node)
        entry[-1] = self.REMOVED
        self.length -= 1
        return entry[0]

    def pop(self):
        """
        Any popped node marked by "REMOVED" does not return, the deleted
        nodes might be popped or still in heap, either case is fine.
        """
        while self.heap:
            priority, node = heapq.heappop(self.heap)
            if node is not self.REMOVED:
                del self.entry_finder[node]
                self.length -= 1
                return priority, node
        raise KeyError('pop from an empty priority queue')
