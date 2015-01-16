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
    for perm in target_perms:
        for block in perm.blocks:
            target_index[block.block_id].append(perm)

    counter = 0
    resolved = _resolve_by_context(ref_perms, repeats)
    to_replace = defaultdict(list)
    for r, contexts in resolved.items():
        one_block_perms = filter(lambda p: len(p.blocks) == 1, target_index[r])
        if len(one_block_perms) != 1: continue

        perm = one_block_perms[0]
        for context in contexts:
            #replace in refs
            to_replace[r].append((context, next_block_id))

            #add new target contig
            new_perm = deepcopy(perm)
            new_perm.blocks[0].block_id = next_block_id
            target_perms.append(new_perm)

            next_block_id += 1
            counter += 1
            #print(r, context)

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
    all_refs = set()

    for perm in ref_perms:
        all_refs.add(perm.genome_name)
        if len(perm.blocks) < 3: continue

        for left_b, mid_b, right_b in zip(perm.blocks[:-2], perm.blocks[1:-1],
                                          perm.blocks[2:]):

            if (mid_b.block_id in repeats and left_b.block_id not in repeats
                                        and right_b.block_id not in repeats):
                left_id = -left_b.signed_id() * mid_b.sign
                right_id = right_b.signed_id() * mid_b.sign
                ref_index[mid_b.block_id][perm.genome_name] \
                                    .append(Context(left_id, right_id))

    resolved_repeats = defaultdict(list)
    for repeat, repeat_index in ref_index.items():
        first_ref = next(iter(all_refs))
        to_compare = repeat_index[first_ref]

        for context in to_compare:
            good = True
            for ref in all_refs:
                if ref == first_ref: continue
                if context not in repeat_index[ref]:
                    good = False
                    break

            if good:
                resolved_repeats[repeat].append(context)

    return resolved_repeats
