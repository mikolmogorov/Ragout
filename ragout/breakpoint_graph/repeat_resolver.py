#(c) 2013-2015 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module resolves repeats so we can
put them into the breakpoint graph.
The underlying intuition is straitforward, however
the code contains a lot of magic :(
"""

from collections import namedtuple, defaultdict
from itertools import chain, product, combinations
from copy import deepcopy, copy
import logging

import networkx as nx

logger = logging.getLogger()


class Context:
    def __init__(self, perm, pos, left, right):
        self.perm = perm
        self.pos = pos
        self.left = left
        self.right = right

    def __str__(self):
        block_id = self.perm.blocks[self.pos].block_id
        return "({0}, {1}, {2}, {3})".format(self.perm.chr_name, block_id,
                                             self.left, self.right)

    def equal(self, other):
        return self.right == other.right and self.left == other.left


def resolve_repeats(ref_perms, target_perms, repeats, phylogeny):
    """
    Does the job
    """
    logger.info("Resolving repeats")

    next_block_id = 0
    for perm in chain(ref_perms, target_perms):
        next_block_id = max(next_block_id,
                            max(map(lambda b: b.block_id, perm.blocks)) + 1)
    first_block_id = next_block_id
    target_name = target_perms[0].genome_name

    ref_contexts = _get_contexts(ref_perms, repeats)
    trg_contexts = _get_contexts(target_perms, repeats)

    #getting matches
    repetitive_matches = []
    unique_matches = []
    for repeat_id, contexts in ref_contexts.items():
        by_genome = defaultdict(list)
        for ctx in contexts:
            by_genome[ctx.perm.genome_name].append(ctx)

        logger.debug("==Resolving {0}".format(repeat_id))
        profiles = _split_into_profiles(by_genome, repeats, phylogeny)
        profiles = list(filter(lambda p: _parsimony_test(p, phylogeny,
                                        target_name), profiles))
        unique_m, repetitive_m = _match_target_contexts(profiles,
                                            trg_contexts[repeat_id], repeats)
        unique_matches.extend(unique_m)
        repetitive_matches.extend(repetitive_m)

        for matches in [unique_m, repetitive_m]:
            for trg_ctx, profile in matches:
                logger.debug("T: {0}".format(trg_ctx))
                for ref_ctx in profile:
                    logger.debug("R: {0}".format(ref_ctx))
                logger.debug("--")
            logger.debug("~~~~~~~~~~~~~~~~")
    ##

    ##resolving unique
    for trg_ctx, profile in unique_matches:
        for ref_ctx in profile:
            assert (trg_ctx.perm.blocks[trg_ctx.pos].block_id ==
                    ref_ctx.perm.blocks[ref_ctx.pos].block_id)
            ref_ctx.perm.blocks[ref_ctx.pos].block_id = next_block_id
        trg_ctx.perm.blocks[trg_ctx.pos].block_id = next_block_id
        next_block_id += 1
    ##

    #resolving repetitive
    by_target_perm = defaultdict(list)
    for trg_ctx, profile in repetitive_matches:
        by_target_perm[trg_ctx.perm].append((trg_ctx, profile))
    to_remove = set()
    for perm, matches in by_target_perm.items():
        groups = _split_by_instance(matches)

        #print(perm)
        for group in groups:
            #for g in group:
                #print(g[0])
                #print(map(str, g[1]))
                #print("--")
            #print(map(lambda p: str(p[0]) + " " + str(p[1]), group))
            new_perm = deepcopy(perm)
            for trg_ctx, profile in group:
                for ref_ctx in profile:
                    assert (new_perm.blocks[trg_ctx.pos].block_id ==
                            ref_ctx.perm.blocks[ref_ctx.pos].block_id)
                    ref_ctx.perm.blocks[ref_ctx.pos].block_id = next_block_id
                new_perm.blocks[trg_ctx.pos].block_id = next_block_id
                next_block_id += 1
            target_perms.append(new_perm)

        if groups:
            to_remove.add(perm)
    target_perms = list(filter(lambda p: p not in to_remove, target_perms))
    ##

    logger.debug("Resolved {0} unique repeat instances"
                        .format(next_block_id - first_block_id))


def _parsimony_test(profile, phylogeny, target_name):
    """
    Determines if the given uniqe instance of a repeat exists in target genome
    """
    states = {g : False for g in phylogeny.terminals_dfs_order()}
    for ctx in profile:
        states[ctx.perm.genome_name] = True

    score_without = phylogeny.estimate_tree(states)
    states[target_name] = True
    score_with = phylogeny.estimate_tree(states)
    return score_with < score_without


def _split_into_profiles(contexts_by_genome, repeats, phylogeny):
    """
    Given repeat contexts in each of reference genomes,
    joins them into "profiles" -- sets of matched contexts
    across different genomes (like an alignemnt column in MSA)
    """
    references = set(contexts_by_genome.keys())
    genomes = filter(lambda g: g in references,
                     phylogeny.terminals_dfs_order())
    profiles  = map(lambda c: [c], contexts_by_genome[genomes[0]])

    #logger.debug(str(genomes))
    for genome in genomes[1:]:
        #finding a matching between existing profiles and a new genome
        genome_ctxs = contexts_by_genome[genome]
        graph = nx.Graph()
        for (pr_id, prof), (ctx_id, ctx) in product(enumerate(profiles),
                                                    enumerate(genome_ctxs)):
            node_prof = "profile" + str(pr_id)
            node_genome = "genome" + str(ctx_id)
            graph.add_node(node_prof, profile=True, prof=prof)
            graph.add_node(node_genome, profile=False, ctx=ctx)

            score = _profile_similarity(prof, ctx, repeats, same_len=True)
            if score > 0:
                graph.add_edge(node_prof, node_genome, weight=score)

        edges = _max_weight_matching(graph)
        for edge in edges:
            prof_node, genome_node = edge
            if graph.node[genome_node]["profile"]:
                prof_node, genome_node = genome_node, prof_node

            #logger.debug("Matched: {0} with score {1}"
            #                .format(graph.node[genome_node]["ctx"],
            #                        graph[prof_node][genome_node]["weight"]))
            #prof, ctx = graph.node[prof_node]["prof"], graph.node[genome_node]["ctx"]
            #for c in prof:
            #    logger.debug("{0}, {1}".format(c, ctx))
            #    logger.debug(_context_similarity(c, ctx, repeats, True))
            #logger.debug(_profile_similarity(graph.node[prof_node]["prof"],
            #                                graph.node[genome_node]["ctx"],
            #                                repeats, True))
            graph.node[prof_node]["prof"].append(graph.node[genome_node]["ctx"])

    return profiles


def _match_target_contexts(profiles, target_contexts, repeats):
    """
    Tries to find a mapping between reference profiles and target contexts
    """
    #TODO: determine if each context exists in target using parsimony procedure
    def is_unique(context):
        return any(b not in repeats for b in
                   map(lambda b: b.block_id, context.perm.blocks))

    unique_matches = []
    repetitive_matches = []

    t_unique = [c for c in target_contexts if is_unique(c)]
    t_repetitive = [c for c in target_contexts if not is_unique(c)]
    #logger.debug("R" + str(map(str, t_repetitive)))

    #create bipartie graph
    graph = nx.Graph()

    #add unique contexts
    for (pr_id, prof), (ctx_id, ctx) in product(enumerate(profiles),
                                                enumerate(t_unique)):
        node_prof = "profile" + str(pr_id)
        node_genome = "target" + str(ctx_id)
        graph.add_node(node_prof, profile=True, prof=prof)
        graph.add_node(node_genome, profile=False, ctx=ctx)

        score = _profile_similarity(prof, ctx, repeats, same_len=False)
        if score > 0:
            graph.add_edge(node_prof, node_genome, weight=score, match="unq")

    #repetetive ones
    different = all(not c_1.equal(c_2) for c_1, c_2 in
                    combinations(t_repetitive, 2))
    if different:
        many_rep = t_repetitive * len(profiles)
        for (pr_id, prof), (ctx_id, ctx) in product(enumerate(profiles),
                                                    enumerate(many_rep)):
            node_prof = "profile" + str(pr_id)
            node_genome = "rep_target" + str(ctx_id)
            graph.add_node(node_prof, profile=True, prof=prof)
            graph.add_node(node_genome, profile=False, ctx=ctx)

            score = _profile_similarity(prof, ctx, repeats, same_len=False)
            if score >= 0:
                graph.add_edge(node_prof, node_genome, weight=score,
                               match="rep")

    edges = _max_weight_matching(graph)
    for edge in edges:
        prof_node, genome_node = edge
        if graph.node[genome_node]["profile"]:
            prof_node, genome_node = genome_node, prof_node

        profile = graph.node[prof_node]["prof"]
        trg_ctx = graph.node[genome_node]["ctx"]

        if graph[prof_node][genome_node]["match"] == "unq":
            unique_matches.append((trg_ctx, profile))
        else:
            repetitive_matches.append((trg_ctx, profile))

        #for ctx in profile:
        #    logger.debug(str(ctx))
        #logger.debug("~~")
        #logger.debug(str(trg_ctx))
        #logger.debug("--")
    #logger.debug("Uniq: {0}, Rep: {1}".format(len(unique_matches),
    #                                          len(repetitive_matches)))

    return unique_matches, repetitive_matches


def _split_by_instance(matches):
    """
    Given matched contexts within a single contig,
    split them into groups where each group corresponds
    to a unique instance of this contig
    """
    by_pos = defaultdict(list)
    trg_ctx_by_pos = {}
    for trg_ctx, profile in matches:
        prof_by_genome = {ctx.perm.genome_name : ctx for ctx in profile}
        by_pos[trg_ctx.pos].append(prof_by_genome)
        trg_ctx_by_pos[trg_ctx.pos] = trg_ctx

    target_perm = matches[0][0].perm
    positions = by_pos.keys()
    #print(target_perm.blocks)

    groups = []
    #try all possible combinations wrt to positions in contig
    for combination in product(*by_pos.values()):
        combination = list(combination)
        master_prof = combination[0]
        master_pos = positions[0]
        try:
            group = [(trg_ctx_by_pos[master_pos], master_prof.values())]
            for genome in master_prof.keys():
                master_ctx = master_prof[genome]
                master_sign = (target_perm.blocks[master_pos].sign *
                               master_ctx.perm.blocks[master_ctx.pos].sign)

                for prof_pos, prof in zip(positions, combination)[1:]:
                    group.append((trg_ctx_by_pos[prof_pos], prof.values()))

                    genome_ctx = prof[genome]
                    if (master_prof[genome].perm.chr_name !=
                        genome_ctx.perm.chr_name):
                        raise KeyError

                    prof_sign = (target_perm.blocks[prof_pos].sign *
                                 genome_ctx.perm.blocks[genome_ctx.pos].sign)
                    shift = (genome_ctx.pos - master_ctx.pos) * prof_sign
                    if (prof_sign != master_sign or
                            shift != prof_pos - master_pos):
                        raise KeyError

            #all good!
            groups.append(group)

        except KeyError:
            continue

    return groups


def _context_similarity(ctx_ref, ctx_trg, repeats, same_len):
    """
    Compute similarity between two contexts
    """
    def alignment(ref, trg):
        """
        Computes global alignment
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
        if same_len:
            for i in xrange(l1):
                table[i][0] = i * GAP
            for i in xrange(l2):
                table[0][i] = i * GAP

        for i, j in product(xrange(1, l1), xrange(1, l2)):
            table[i][j] = max(table[i-1][j] + GAP, table[i][j-1] + GAP,
                              table[i-1][j-1] + match(ref[i-1], trg[j-1]))
        return table[-1][-1]

    if len(ctx_trg.left) + len(ctx_trg.right) == 0:
        return 0

    left = alignment(ctx_ref.left, ctx_trg.left)
    right = alignment(ctx_ref.right[::-1], ctx_trg.right[::-1])
    return left + right


def _profile_similarity(profile, genome_ctx, repeats, same_len):
    """
    Compute similarity of set of contexts vs one context
    """
    scores = list(map(lambda c: _context_similarity(c, genome_ctx,
                                    repeats, same_len), profile))
    return float(sum(scores)) / len(scores)


def _max_weight_matching(graph):
    edges = nx.max_weight_matching(graph, maxcardinality=True)
    unique_edges = set()
    for v1, v2 in edges.items():
        if not (v2, v1) in unique_edges:
            unique_edges.add((v1, v2))

    return list(unique_edges)


def _get_contexts(permutations, repeats):
    """
    Get repeats' contexts
    """
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
