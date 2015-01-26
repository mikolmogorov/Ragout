#(c) 2013-2015 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module tries to resolve simple repeats so we are able to
put them into breakpoint graph
"""

from collections import namedtuple, defaultdict
from itertools import chain, product
from copy import deepcopy, copy
import logging

import networkx as nx

#Context = namedtuple("Context", ["perm", "pos", "left", "right"])
class Context:
    def __init__(self, perm, pos, left, right):
        self.perm = perm
        self.pos = pos
        self.left = left
        self.right = right

    def __str__(self):
        block_id = self.perm.blocks[self.pos].block_id
        return "({0}, {1}, {2})".format(self.perm.chr_name,
                                        self.left, self.right)

logger = logging.getLogger()

def resolve_repeats(ref_perms, target_perms, repeats):
    """
    Does the job
    """
    logger.info("Resolving repeats")

    next_block_id = 0
    for perm in chain(ref_perms, target_perms):
        next_block_id = max(next_block_id,
                            max(map(lambda b: b.block_id, perm.blocks)) + 1)

    ref_contexts = _get_contexts(ref_perms, repeats)
    trg_contexts = _get_contexts(target_perms, repeats)
    one2one, one2many = _match_contexts(ref_contexts, trg_contexts)

    for trg_ctx, ref_ctx in one2one:
        #logger.debug("{0} with {1}".format(trg_ctx, ref_ctx))
        trg_ctx.perm.blocks[trg_ctx.pos].block_id = next_block_id
        ref_ctx.perm.blocks[ref_ctx.pos].block_id = next_block_id
        next_block_id += 1


    for trg_perm, ref_ctx in one2many:
        #logger.debug(str(ref_ctx))
        assert len(trg_perm.blocks) == 1
        new_perm = deepcopy(trg_perm)
        new_perm.blocks[0].block_id = next_block_id
        target_perms.append(new_perm)
        ref_ctx.perm.blocks[ref_ctx.pos].block_id = next_block_id
        next_block_id += 1


def _edit_distance(ref, trg):
    """
    Computes edit distance, allowing
    free gaps on the left side
    """
    GAP = 1
    MISS = 1
    l1, l2 = len(ref) + 1, len(trg) + 1
    table = [[0 for _ in xrange(l2)] for _ in xrange(l1)]

    """
    for i in xrange(l1):
        table[i][0] = i
    for j in xrange(l2):
        table[0][j] = j
    """

    for i, j in product(xrange(1, l1), xrange(1, l2)):
        table[i][j] = min(table[i-1][j] + GAP, table[i][j-1] + GAP,
                          table[i-1][j-1] + (MISS if ref[i-1] != trg[j-1]
                                                  else 0))
    return table[-1][-1]


def _context_distance(ctx_ref, ctx_trg):
    """
    Compute discordance between two contexts
    """
    length = len(ctx_trg.left) + len(ctx_trg.right)
    if not length: return 0

    left = _edit_distance(ctx_ref.left, ctx_trg.left)
    right = _edit_distance(ctx_ref.right[::-1], ctx_trg.right[::-1])
    return float(left + right) / length


def _min_weight_matching(graph):
    for v1, v2 in graph.edges_iter():
        graph[v1][v2]["weight"] = -graph[v1][v2]["weight"] #want minimum weight

    edges = nx.max_weight_matching(graph, maxcardinality=True)
    unique_edges = set()
    for v1, v2 in edges.items():
        if not (v2, v1) in unique_edges:
            unique_edges.add((v1, v2))

    return list(unique_edges)


def _match_contexts(ref_contexts, target_contexts):
    MAX_DIST = 0.4
    matched_ref_ctx = ref_contexts      #now let's assume we have one reference

    one2one_matched = []
    one2many_matched = []

    for block in target_contexts:
        t_contexts = target_contexts[block]
        r_contexts = ref_contexts[block]

        t_partial = [c for c in t_contexts if len(c.left) + len(c.right) > 0]
        t_zero = [c for c in t_contexts if len(c.left) + len(c.right) == 0]

        #create bipartie graph
        logger.debug("Processing {0}".format(block))
        logger.debug("Target contexts:\n{0}"
                            .format("\n".join(map(str, t_contexts))))
        logger.debug("Reference contexts:\n{0}"
                            .format("\n".join(map(str, r_contexts))))
        graph = nx.Graph()
        for (no_t, ctx_t), (no_r, ctx_r) in product(enumerate(t_partial),
                                                    enumerate(r_contexts)):
            node_ref, node_trg = "ref" + str(no_r), "trg" + str(no_t)
            graph.add_node(node_ref, ref=True, ctx=ctx_r)
            graph.add_node(node_trg, ref=False, ctx=ctx_t)

            distance = _context_distance(ctx_r, ctx_t)
            #logger.debug("D {0} -- {1} {2}".format(ctx_r, ctx_t, distance))
            if distance < MAX_DIST:
                graph.add_edge(node_ref, node_trg, weight=distance)

        edges = _min_weight_matching(graph)
        used_contexts = set()
        for edge in edges:
            u, v = edge
            if graph.node[u]["ref"]:
                u, v = v, u
            one2one_matched.append((graph.node[u]["ctx"],
                                    graph.node[v]["ctx"]))
            used_contexts.add(graph.node[v]["ctx"])
            logger.debug("M: {0} -- {1}".format(graph.node[v]["ctx"],
                                                graph.node[u]["ctx"]))

        if len(t_zero) == 1:
            #for node in graph.nodes():
            for r_ctx in r_contexts:
                if r_ctx not in used_contexts:
                    one2many_matched.append((t_zero[0].perm, r_ctx))
                    logger.debug("Z {0} -- {1}".format(t_zero[0].perm.chr_name,
                                                       r_ctx))
                    #print("ya")

    return one2one_matched, one2many_matched


"""
def _match_contexts(ref_contexts, target_contexts):
    MIN_DIST = 0.4
    EPS = 0.01

    matched_ref_ctx = ref_contexts      #now let's assume we have one reference

    one2one_matches = []
    zero_ctx_matches = []

    for block, contexts in target_contexts.items():
        partial_ctx = [c for c in contexts if len(c.left) + len(c.right) > 0]
        zero_ctx = [c for c in contexts if len(c.left) + len(c.right) == 0]

        #contexts_left = matched_ref_ctx[block][::]
        partial_left = False
        for trg_ctx in partial_ctx:
            candidates = []
            for ref_ctx in matched_ref_ctx[block]:
                distance = _context_distance(ref_ctx, trg_ctx)
                if distance < MIN_DIST:
                    candidates.append((ref_ctx, distance))

            if not candidates:
                partial_left = True
                continue

            candidates.sort(key=lambda c: c[1])
                #logger.debug("Resolving {0}".format(block))
                #logger.debug("No match for {0}".format(trg_ctx))
            if len(candidates) == 1:
                one2one_matches.append((candidates[0][0], trg_ctx))
                contexts_left.remove(candidates[0][0])
                #logger.debug("{0} matched with {1}"
                #                    .format(trg_ctx, candidates[0]))
            if len(candidates) > 1:
                if abs(candidates[0][1] - candidates[1][1]) > EPS:
                    one2one_matches.append((candidates[0][0], trg_ctx))
                    #contexts_left.remove(candidates[0][0])
                #else:
                #    for cand in candidates:
                #        contexts_left.remove(cand[0])

        if len(zero_ctx) == 1 and not partial_left:
            #print(len(contexts_left))
            print("let's have some fun!")
                #else:
                    #logger.debug("Resolving {0}".format(block))
                    #logger.debug("{0} has multiple candidates: {1}"
                    #                    .format(trg_ctx, candidates))

    #logger.debug("Resolved {0} contexts".format(resolved_count))
    #logger.debug("Unresolved {0} contexts".format(unresolved_count))
"""


def _get_contexts(permutations, repeats):
    WINDOW = 5

    contexts = defaultdict(list)
    for perm in permutations:
        for pos in xrange(len(perm.blocks)):
            block = perm.blocks[pos]
            if block.block_id not in repeats:
                continue

            left_start = max(0, pos - WINDOW)
            left_end = max(0, pos)
            left_context = list(map(lambda b: b.signed_id() * block.sign,
                                    perm.blocks[left_start:left_end]))

            right_start = min(len(perm.blocks), pos + 1)
            right_end = min(len(perm.blocks), pos + WINDOW + 1)
            right_context = list(map(lambda b: b.signed_id() * block.sign,
                                     perm.blocks[right_start:right_end]))

            if block.sign < 0:
                left_context, right_context = (right_context[::-1],
                                               left_context[::-1])

            contexts[block.block_id].append(Context(perm, pos, left_context,
                                                    right_context))
    return contexts
