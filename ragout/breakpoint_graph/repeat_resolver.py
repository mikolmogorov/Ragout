#(c) 2013-2015 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

from collections import namedtuple, defaultdict

def resolve_repeats(ref_perms, target_perms, repeats):
    good_repeats = find_good(ref_perms, repeats)
    target_index = defaultdict(list)
    for perm in target_perms:
        for block in perm.blocks:
            target_index[block.block_id].append(perm)

    for r, context in good_repeats.items():
        for perm in target_index[r]:
            if len(perm.blocks) == 1:
                print(r, context)
                print(list(map(lambda b: b.signed_id(), perm.blocks)))

def find_good(ref_perms, repeats):
    ref_index = defaultdict(lambda: defaultdict(list))
    all_refs = set()

    for perm in ref_perms:
        all_refs.add(perm.genome_name)
        if len(perm.blocks) < 3:
            continue

        for left_b, mid_b, right_b in zip(perm.blocks[:-2], perm.blocks[1:-1],
                                          perm.blocks[2:]):

            if (mid_b.block_id in repeats and left_b.block_id not in repeats
                                        and right_b.block_id not in repeats):
                left_id = -left_b.signed_id() * mid_b.sign
                right_id = right_b.signed_id() * mid_b.sign
                ref_index[mid_b.block_id][perm.genome_name].append((left_id,
                                                                    right_id))

    good_repeats = defaultdict(list)
    #now find some good guys!
    for repeat, repeat_index in ref_index.items():
        first_ref = next(iter(all_refs))
        to_compare = repeat_index[first_ref]

        for context in to_compare:
            good = True
            for ref in all_refs:
                if ref == first_ref:
                    continue
                if context not in repeat_index[ref]:
                    good = False
                    break

            if good:
                good_repeats[repeat].append(context)

    return good_repeats
