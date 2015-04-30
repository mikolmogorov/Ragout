#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module merges scaffolds from two consecutive
iterations
"""

from collections import namedtuple, defaultdict
from itertools import product, chain
import sys
import logging

from ragout.shared.datatypes import Contig, Scaffold

logger = logging.getLogger()


def merge(big_scaffolds, small_scaffolds):
    """
    Merges two assemblies assuming that big one is correct
    """
    count_diff_scaf = 0
    count_diff_orient = 0
    count_inconsistent = 0

    total_success = 0
    total_fail = 0
    total_inserted = 0
    not_found = 0

    logger.info("Merging two iterations")
    big_count = defaultdict(int)
    for scf in big_scaffolds:
        for c in scf.contigs:
            big_count[c.seq_name] += 1

    small_count = defaultdict(int)
    for scf in small_scaffolds:
        for c in scf.contigs:
            small_count[c.seq_name] += 1

    repeats = set(seq for (seq, count) in
                  chain(big_count.items(), small_count.items()) if count > 1)
    big_unique = set(seq for (seq, count) in big_count.items() if count == 1)

    small_index = {}
    for scf in small_scaffolds:
        for pos, contig in enumerate(scf.contigs):
            if contig.seq_name not in repeats:
                assert contig.seq_name not in small_index
                small_index[contig.seq_name] = (scf, pos)

    new_scafflods = []
    for big_scf in big_scaffolds:
        new_contigs = []
        non_repeats = list(filter(lambda i: big_scf.contigs[i].seq_name
                                        not in repeats,
                                  xrange(len(big_scf.contigs))))
        for left_idx, right_idx in zip(non_repeats[:-1], non_repeats[1:]):
            left_cnt = big_scf.contigs[left_idx]
            right_cnt = big_scf.contigs[right_idx]

            consistent = False
            if (left_cnt.seq_name in small_index and
                right_cnt.seq_name in small_index):
                consistent = True
                left_scf, left_pos = small_index[left_cnt.seq_name]
                right_scf, right_pos = small_index[right_cnt.seq_name]

                big_sign = left_cnt.sign == right_cnt.sign
                small_sign = (left_scf.contigs[left_pos].sign ==
                              right_scf.contigs[right_pos].sign)

                if left_scf != right_scf:
                    count_diff_scaf += 1
                    consistent = False
                elif big_sign != small_sign:
                    count_diff_orient += 1
                    consistent = False
                else:
                    same_dir = left_pos < right_pos
                    if not same_dir:
                        left_pos, right_pos = right_pos, left_pos

                    weak_contigs = left_scf.contigs[left_pos + 1 : right_pos]
                    if any(c.seq_name in big_unique for c in weak_contigs):
                        count_inconsistent += 1
                        consistent = False

                    if not same_dir:
                        weak_contigs = list(map(lambda c: c.reverse_copy(),
                                                weak_contigs[::-1]))
                    link_to_change = left_scf.contigs[left_pos].link
            else:
                not_found += 1

            new_contigs.append(left_cnt)
            if consistent:
                new_contigs[-1].link = link_to_change
                new_contigs.extend(weak_contigs)
                total_success += 1
                total_inserted += len(weak_contigs)
                #logger.debug("Inserting '{0}' between {1} and {2}"
                #             .format(map(lambda c: c.seq_name, weak_contigs),
                #                     left_cnt, right_cnt))
            else:
                new_contigs.extend(big_scf.contigs[left_idx+1:right_idx])
                total_fail += 1

        if len(new_contigs) > 1:
            new_contigs.append(right_cnt)
            s = Scaffold(big_scf.name)
            s.contigs = new_contigs
            new_scafflods.append(s)
        else:   #because of repeats
            new_scafflods.append(big_scf)

    logger.debug("Fail: not found: {0}".format(not_found))
    logger.debug("Fail: different scaffolds: {0}".format(count_diff_scaf))
    logger.debug("Fail: different orientatilns: {0}".format(count_diff_orient))
    logger.debug("Fail: inconsistent: {0}".format(count_inconsistent))
    logger.debug("Total success: {0}".format(total_success))
    logger.debug("Total fail: {0}".format(total_fail))
    logger.debug("Total inserted: {0}".format(total_inserted))

    num_contigs = 0
    for scf in new_scafflods:
        num_contigs += len(scf.contigs)
    logger.debug("Result: {0} contigs in {1} scaffolds"
                                    .format(num_contigs, len(new_scafflods)))

    return new_scafflods
