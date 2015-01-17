#(c) 2013-2015 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module tries to resolve simple repeats so we are able to
put them into breakpoint graph
"""

from collections import namedtuple, defaultdict
from itertools import chain
from copy import deepcopy
import logging

Context = namedtuple("Context", ["left", "right"])
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

    target_index = defaultdict(list)
    target_adj = set()
    for perm in target_perms:
        for block in perm.blocks:
            target_index[block.block_id].append(perm)
        if len(perm.blocks) > 1:
            for lb, rb in zip(perm.blocks[:-1], perm.blocks[1:]):
                target_adj.add(-lb.signed_id())
                target_adj.add(rb.signed_id())

    counter = 0
    resolved = _resolve_by_context(ref_perms, repeats)
    to_replace = defaultdict(list)
    for r, contexts in resolved.items():
        if len(target_index[r]) != 1: continue
        perm = target_index[r][0]
        if len(perm.blocks) != 1: continue

        for context in contexts:
            if len(set(context).intersection(target_adj)): continue

            #replace in refs
            to_replace[r].append((context, next_block_id))

            #add new target contig
            new_perm = deepcopy(perm)
            new_perm.blocks[0].block_id = next_block_id
            target_perms.append(new_perm)

            next_block_id += 1
            counter += 1

    _replace_blocks(to_replace, ref_perms)
    logger.debug("{0} repeat instances resolved".format(counter))


def _replace_blocks(to_replace, permutations):
    """
    Replaces block instance with a given context in all permutations
    """
    for perm in permutations:
        for left_b, mid_b, right_b in zip(perm.blocks[:-2], perm.blocks[1:-1],
                                          perm.blocks[2:]):

            left_id = -left_b.signed_id() * mid_b.sign
            right_id = right_b.signed_id() * mid_b.sign
            for context, new_id in to_replace[mid_b.block_id]:
                if context == (left_id, right_id):
                    mid_b.block_id = new_id


def _resolve_by_context(ref_perms, repeats):
    """
    Tries to resolve trivial repeats (that have same context in every ref)
    """
    ref_index = defaultdict(lambda: defaultdict(list))
    ref_with_block = defaultdict(set)
    for perm in ref_perms:
        for block in perm.blocks:
            ref_with_block[block.block_id].add(perm.genome_name)

    #collecting all repeat contexts
    for perm in ref_perms:
        if len(perm.blocks) < 3: continue

        for left_b, mid_b, right_b in zip(perm.blocks[:-2], perm.blocks[1:-1],
                                          perm.blocks[2:]):

            if (mid_b.block_id in repeats and left_b.block_id not in repeats
                and right_b.block_id not in repeats):
                left_id = -left_b.signed_id() * mid_b.sign
                right_id = right_b.signed_id() * mid_b.sign
                ref_index[mid_b.block_id][perm.genome_name] \
                                    .append(Context(left_id, right_id))

    resolved_repeats = defaultdict(set)
    for repeat, repeat_index in ref_index.items():
        refs_to_check = ref_with_block[repeat]
        #first ref and context
        for ref in refs_to_check:
            for context in repeat_index[ref]:
                other_refs = (refs_to_check
                             .intersection(ref_with_block[abs(context.left)])
                             .intersection(ref_with_block[abs(context.right)]))

                good = True
                for other_ref in other_refs:
                    if other_ref == ref: continue
                    if context not in repeat_index[other_ref]:
                        good = False
                        break

                if good:
                    resolved_repeats[repeat].add(context)

    return resolved_repeats
