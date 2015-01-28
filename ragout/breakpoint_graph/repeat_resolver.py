#(c) 2013-2015 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module tries to resolve repeats so we are able to
put them into breakpoint graph
"""

from collections import namedtuple, defaultdict
from itertools import chain, product, combinations
from copy import deepcopy, copy
import logging

import networkx as nx

class Context:
    def __init__(self, perm, pos, left, right):
        self.perm = perm
        self.pos = pos
        self.left = left
        self.right = right

    def __str__(self):
        block_id = self.perm.blocks[self.pos].block_id
        return "({0}, {1}, {2}, {3})".format(self.perm.chr_name, self.pos,
                                             self.left, self.right)

    def equal(self, other):
        return self.right == other.right and self.left == other.left


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
    unique_matches, repetitive_matches = _match_contexts(ref_contexts,
                                                         trg_contexts, repeats)

    for trg_ctx, ref_ctx in unique_matches:
        assert (trg_ctx.perm.blocks[trg_ctx.pos].block_id ==
                ref_ctx.perm.blocks[ref_ctx.pos].block_id)
        trg_ctx.perm.blocks[trg_ctx.pos].block_id = next_block_id
        ref_ctx.perm.blocks[ref_ctx.pos].block_id = next_block_id
        next_block_id += 1

    #now processing repetitive contigs
    by_target = defaultdict(list)
    for trg_ctx, ref_ctx in repetitive_matches:
        by_target[trg_ctx.perm].append((trg_ctx, ref_ctx))

    for perm, contexts in by_target.items():
        if any(b.block_id not in repeats for b in perm.blocks):
            continue

        logger.debug("Perm: {0}".format(perm.chr_name))
        groups = _split_by_instance(contexts)

        for group in groups:
            new_perm = deepcopy(perm)
            for trg_ctx, ref_ctx in group:
                assert (new_perm.blocks[trg_ctx.pos].block_id ==
                        ref_ctx.perm.blocks[ref_ctx.pos].block_id)
                new_perm.blocks[trg_ctx.pos].block_id = next_block_id
                ref_ctx.perm.blocks[ref_ctx.pos].block_id = next_block_id
                next_block_id += 1

            target_perms.append(new_perm)

        target_perms.remove(perm)


def _split_by_instance(contexts):
    """
    Given a set of matched contexts, group them by instances of
    coresponding contig in target genome
    """
    index = defaultdict(lambda: defaultdict(list))
    for trg_ctx, ref_ctx in contexts:
        index[ref_ctx.perm][trg_ctx.pos].append((trg_ctx, ref_ctx))

    groups = []
    for ref_perm, by_pos in index.items():
        if len(by_pos) != len(contexts[0][0].perm.blocks): continue

        logger.debug("Ref: {0}".format(ref_perm.chr_name))

        #some magic position shifting
        positions = {}
        zero_pos = sorted(by_pos.keys())[0]
        for pos, contexts in by_pos.items():
            shift = pos - zero_pos
            sign = lambda t, r: (t.perm.blocks[t.pos].sign *
                                 r.perm.blocks[r.pos].sign)
            positions[pos] = list(map(lambda (t, r): r.pos - shift * sign(t, r),
                                      contexts))
            logger.debug("Pos {0} : {1}".format(pos, positions[pos]))

        #now find equal positions
        for num, ref_pos in enumerate(positions[zero_pos]):
            contexts_group = [by_pos[zero_pos][num]]
            try:
                for pos in positions:
                    if pos == zero_pos: continue
                    index = positions[pos].index(ref_pos)
                    contexts_group.append(by_pos[pos][index])
                groups.append(contexts_group)
                logger.debug("Group:\n{0}"
                    .format("\n".join(map(lambda (t, r): str((str(t), str(r))),
                                      contexts_group))))
            except ValueError:
                pass

    return groups


def _alignment(ref, trg, repeats):
    """
    Computes global alignment, allowing
    free gaps on the left side
    """
    GAP = -1
    def match(a, b):
        if a != b:
            return -2
        if abs(a) in repeats:
            return 1
        return 2

    l1, l2 = len(ref) + 1, len(trg) + 1
    table = [[0 for _ in xrange(l2)] for _ in xrange(l1)]

    for i, j in product(xrange(1, l1), xrange(1, l2)):
        table[i][j] = max(table[i-1][j] + GAP, table[i][j-1] + GAP,
                          table[i-1][j-1] + match(ref[i-1], trg[j-1]))
    return table[-1][-1]


def _similarity_score(ctx_ref, ctx_trg, repeats):
    """
    Compute similarity between two contexts
    """
    if len(ctx_trg.left) + len(ctx_trg.right) == 0:
        return 0

    left = _alignment(ctx_ref.left, ctx_trg.left, repeats)
    right = _alignment(ctx_ref.right[::-1], ctx_trg.right[::-1], repeats)
    return left + right


def _max_weight_matching(graph):
    edges = nx.max_weight_matching(graph, maxcardinality=True)
    unique_edges = set()
    for v1, v2 in edges.items():
        if not (v2, v1) in unique_edges:
            unique_edges.add((v1, v2))

    return list(unique_edges)


def _match_contexts(ref_contexts, target_contexts, repeats):
    def is_unique(context):
        return any(abs(b) not in repeats for b in
                              chain(context.left, context.right))

    unique_matches = []
    repetitive_matches = []

    for block in target_contexts:
        t_contexts = target_contexts[block]
        r_contexts = ref_contexts[block]

        t_unique = [c for c in t_contexts if is_unique(c)]
        t_repetitive = [c for c in t_contexts if not is_unique(c)]

        #create bipartie graph
        logger.debug("Processing {0}".format(block))
        logger.debug("Target contexts:\n{0}"
                            .format("\n".join(map(str, t_contexts))))
        logger.debug("Reference contexts:\n{0}"
                            .format("\n".join(map(str, r_contexts))))
        graph = nx.Graph()

        #add unique contexts
        for (no_t, ctx_t), (no_r, ctx_r) in product(enumerate(t_unique),
                                                    enumerate(r_contexts)):
            node_ref, node_trg = "ref" + str(no_r), "unq" + str(no_t)
            graph.add_node(node_ref, ref=True, ctx=ctx_r)
            graph.add_node(node_trg, ref=False, ctx=ctx_t)

            score = _similarity_score(ctx_r, ctx_t, repeats)
            if score > 0:
                graph.add_edge(node_ref, node_trg, weight=score, match="unq")

        #repetetive ones
        different = all(not c_1.equal(c_2) for c_1, c_2 in
                                           combinations(t_repetitive, 2))
        if different:
            logger.debug("Repetetive: {0}".format(map(str, t_repetitive)))
            many_rep = t_repetitive * len(r_contexts)
            for (no_t, ctx_t), (no_r, ctx_r) in product(enumerate(many_rep),
                                                        enumerate(r_contexts)):
                node_ref, node_trg = "ref" + str(no_r), "rep" + str(no_t)
                graph.add_node(node_ref, ref=True, ctx=ctx_r)
                graph.add_node(node_trg, ref=False, ctx=ctx_t)

                score = _similarity_score(ctx_r, ctx_t, repeats)
                if score >= 0:
                    graph.add_edge(node_ref, node_trg, weight=score,
                                   match="rep")
        else:
            logger.debug("Repetitive contexts are ambigous")

        edges = _max_weight_matching(graph)
        for edge in edges:
            trg_node, ref_node = edge
            if graph.node[trg_node]["ref"]:
                trg_node, ref_node = ref_node, trg_node

            match_pair = (graph.node[trg_node]["ctx"],
                          graph.node[ref_node]["ctx"])
            if graph[trg_node][ref_node]["match"] == "unq":
                unique_matches.append(match_pair)
                logger.debug("M: {0} -- {1}".format(graph.node[trg_node]["ctx"],
                                                   graph.node[ref_node]["ctx"]))
            else:
                repetitive_matches.append(match_pair)
                logger.debug("Z: {0} -- {1}".format(graph.node[ref_node]["ctx"],
                                                   graph.node[ref_node]["ctx"]))

    return unique_matches, repetitive_matches


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
