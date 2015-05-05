#(c) 2013-2015 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module resolves repeats so we can
put them into the breakpoint graph.
The underlying intuition is straitforward, however
the code contains a lot of magic. Sorry :(
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
        return "({0}, {1}, {2}, {3})".format(self.perm.chr_name, self.pos,
                                             self.left, self.right)

    def equal(self, other):
        return self.right == other.right and self.left == other.left

MP = namedtuple("MatchPair", ["trg", "prof"])
class MatchPair(MP):
    def __hash__(self):
        return id(self)


def resolve_repeats(ref_perms, target_perms, repeats,
                    phylogeny, draft_refs):
    """
    Does the job
    """
    logger.info("Resolving repeats")
    logger.debug("Unique repeat blocks: {0}".format(len(repeats)))

    next_block_id = 0
    for perm in chain(ref_perms, target_perms):
        next_block_id = max(next_block_id,
                            max(map(lambda b: b.block_id, perm.blocks)) + 1)
    first_block_id = next_block_id
    target_name = target_perms[0].genome_name

    ref_contexts = _get_contexts(ref_perms, repeats)
    trg_contexts = _get_contexts(target_perms, repeats)

    purely_repetitive = 0
    for perm in target_perms:
        if all(map(lambda b: b.block_id in repeats, perm.blocks)):
            purely_repetitive += 1
    logger.debug("Purely repetitive sequences: {0}".format(purely_repetitive))

    #getting matches
    repetitive_matches = []
    unique_matches = []
    for repeat_id, contexts in ref_contexts.items():
        by_genome = defaultdict(list)
        for ctx in contexts:
            by_genome[ctx.perm.genome_name].append(ctx)

        #logger.debug("==Resolving {0}".format(repeat_id))
        profiles = _split_into_profiles(by_genome, repeats, phylogeny)
        parsimony_test = lambda p: _parsimony_test(p, phylogeny, target_name,
                                                   draft_refs)
        profiles = list(filter(parsimony_test, profiles))
        unique_m, repetitive_m = _match_target_contexts(profiles,
                                            trg_contexts[repeat_id], repeats)
        unique_matches.extend(unique_m)
        repetitive_matches.extend(repetitive_m)

        #for matches in [unique_m, repetitive_m]:
        #    for trg_ctx, profile in matches:
        #        logger.debug("T: {0}".format(trg_ctx))
        #        for ref_ctx in profile:
        #            logger.debug("R: {0}".format(ref_ctx))
        #        logger.debug("--")
        #    logger.debug("~~~~~~~~~~~~~~~~")
    ##

    matched_contigs = set()
    for m in repetitive_matches:
        if all(map(lambda b: b.block_id in repeats, m.trg.perm.blocks)):
            matched_contigs.add(m.trg.perm)
    logger.debug("Repetitive sequences with matches: {0}".format(len(matched_contigs)))

    #logger.debug("Discarded: {0}".format(g_discarded))
    #logger.debug("Unmatched: {0}".format(g_unmatched))
    logger.debug("Unique matches: {0}".format(len(unique_matches)))
    logger.debug("Repetitive matches: {0}".format(len(repetitive_matches)))

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
    for match in repetitive_matches:
        by_target_perm[match.trg.perm].append(match)
    to_remove = set()
    new_contigs = 0
    for perm, matches in by_target_perm.items():
        groups = _split_by_instance(matches)

        for group in groups:
            new_perm = deepcopy(perm)
            for trg_ctx, profile in group:
                for ref_ctx in profile:
                    assert (new_perm.blocks[trg_ctx.pos].block_id ==
                            ref_ctx.perm.blocks[ref_ctx.pos].block_id)
                    ref_ctx.perm.blocks[ref_ctx.pos].block_id = next_block_id
                new_perm.blocks[trg_ctx.pos].block_id = next_block_id
                next_block_id += 1
            target_perms.append(new_perm)
            new_contigs += 1

        if groups:
            to_remove.add(perm)
    target_perms = list(filter(lambda p: p not in to_remove, target_perms))
    ##

    logger.debug("Resolved {0} unique repeat instances"
                        .format(next_block_id - first_block_id))
    logger.debug("Saved sequences: {0}".format(len(to_remove)))
    logger.debug("Added {0} extra contigs".format(new_contigs))


def _parsimony_test(profile, phylogeny, target_name, draft_refs):
    """
    Determines if the given uniqe instance of a repeat exists in target genome
    """
    states = {g : False if g not in draft_refs else None
              for g in phylogeny.terminals_dfs_order()}
    for ctx in profile:
        states[ctx.perm.genome_name] = True

    states[target_name] = False
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

            graph.node[prof_node]["prof"].append(graph.node[genome_node]["ctx"])

    return profiles

def _match_target_contexts(profiles, target_contexts, repeats):
    """
    Tries to find a mapping between reference profiles and target contexts
    """
    def is_unique(context):
        return any(b not in repeats for b in
                   map(lambda b: b.block_id, context.perm.blocks))

    unique_matches = []
    repetitive_matches = []

    t_unique = [c for c in target_contexts if is_unique(c)]
    t_repetitive = [c for c in target_contexts if not is_unique(c)]

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
    dups = set()
    for ctx_1, ctx_2 in combinations(t_repetitive, 2):
        if ctx_1.equal(ctx_2):
            dups.add(ctx_1)
            dups.add(ctx_2)
    different = [ctx for ctx in t_repetitive if ctx not in dups]
    #if dups:
    #    logger.debug(map(str, dups))

    if different:
        many_rep = different * len(profiles)
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

    #get matching
    target_matched = set()
    edges = _max_weight_matching(graph)
    for edge in edges:
        prof_node, genome_node = edge
        if graph.node[genome_node]["profile"]:
            prof_node, genome_node = genome_node, prof_node

        profile = graph.node[prof_node]["prof"]
        trg_ctx = graph.node[genome_node]["ctx"]

        if graph[prof_node][genome_node]["match"] == "unq":
            unique_matches.append(MatchPair(trg_ctx, profile))
        else:
            repetitive_matches.append(MatchPair(trg_ctx, profile))
            target_matched.add(trg_ctx)

    """
    if different:
        g_unmatched += len(t_repetitive) - len(target_matched)
        for cxt in t_repetitive:
            if cxt not in target_matched:
                logger.debug("{0}".format(str(cxt)))
        if len(t_repetitive) - len(target_matched) > 0:
            logger.debug("trg: {0}, prof: {1}".format(len(target_contexts),
                                                      len(profiles)))
            logger.debug(map(str, target_contexts))
    else:
        g_discarded += len(t_repetitive)
    """

    return unique_matches, repetitive_matches


def _split_by_instance(matches):
    """
    Given matched contexts within a single contig,
    split them into groups where each group corresponds
    to a unique instance of this contig
    """
    target_perm = matches[0][0].perm
    if len(target_perm.blocks) == 1:    #trivial case
        return list(map(lambda m: [m], matches))

    #logger.debug("=========")
    #logger.debug(target_perm)
    #logger.debug("=========")

    by_pos = defaultdict(list)
    for match in matches:
        by_pos[match.trg.pos].append(match)
    positions = sorted(by_pos.keys())

    def prof_agreement(match_1, match_2):
        index_1 = {ctx.perm.genome_name : ctx for ctx in match_1.prof}
        index_2 = {ctx.perm.genome_name : ctx for ctx in match_2.prof}
        shared_genomes = set(index_1.keys()) & set(index_2.keys())
        if not len(shared_genomes):
            return False
        for genome in shared_genomes:
            if index_1[genome].perm.chr_name != index_2[genome].perm.chr_name:
                return False
            sign_1 = (target_perm.blocks[match_1.trg.pos].sign *
                      index_1[genome].perm.blocks[index_1[genome].pos].sign)
            sign_2 = (target_perm.blocks[match_2.trg.pos].sign *
                      index_2[genome].perm.blocks[index_2[genome].pos].sign)
            if sign_1 != sign_2:
                return False
            shift = (index_2[genome].pos - index_1[genome].pos) * sign_1
            if shift != match_2.trg.pos - match_1.trg.pos:
                return False
        return True

    groups = map(lambda x: [x], by_pos[positions[0]])
    #now try to extend each group
    for pos in positions[1:]:
        unused_matches = set(by_pos[pos])
        for group in groups:
            prev_match = group[-1]
            for next_match in by_pos[pos]:
                if (next_match in unused_matches and
                        prof_agreement(prev_match, next_match)):
                    unused_matches.remove(next_match)
                    group.append(next_match)
                    break
        groups.extend([[m] for m in unused_matches])

    #min_group = max(2, (len(target_perm.blocks) + 1) / 2)
    min_group = len(target_perm.blocks) / 2 + 1
    #for group in groups:
    #    if len(group) < min_group:
    #        continue
    #    logger.debug("##group")
    #    for match in group:
    #        logger.debug("####profile pos {0}".format(match.trg.pos))
    #        for ctx in match.prof:
    #            logger.debug(ctx)
    #    logger.debug("")

    return list(filter(lambda g: len(g) >= min_group, groups))


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
