#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module merges scaffolds from two consecutive
iterations
"""

from collections import namedtuple, defaultdict
from itertools import product
import sys
import logging

from ragout.shared.datatypes import Contig, Scaffold

logger = logging.getLogger()


def merge(big_scaffolds, small_scaffolds):
    """
    The only function here
    """
    logger.info("Merging two iterations")
    big_index = set()
    for scf in big_scaffolds:
        for c in scf.contigs:
            big_index.add(c.seq_name)

    small_index = defaultdict(list)
    for scf in small_scaffolds:
        for pos, contig in enumerate(scf.contigs):
            small_index[contig.seq_name].append((scf, pos))

    new_scafflods = []
    for scf in big_scaffolds:
        new_contigs = []
        for left_cnt, right_cnt in zip(scf.contigs[:-1], scf.contigs[1:]):
            new_contigs.append(left_cnt)

            left_candidates = small_index[left_cnt.seq_name]
            right_candidates = small_index[right_cnt.seq_name]

            #handling repeats too
            consistent = []
            for left_cand, right_cand in product(left_candidates,
                                                 right_candidates):
                left_scf, left_pos = left_cand
                right_scf, right_pos = right_cand

                big_sign = left_cnt.sign == right_cnt.sign
                small_sign = (left_scf.contigs[left_pos].sign ==
                              right_scf.contigs[right_pos].sign)

                if (left_scf != right_scf or
                        abs(left_pos - right_pos) == 1 or
                        big_sign != small_sign):
                    continue

                same_dir = left_pos < right_pos
                if not same_dir:
                    left_pos, right_pos = right_pos, left_pos

                weak_contigs = left_scf.contigs[left_pos + 1 : right_pos]
                if any(c in big_index for c in weak_contigs):
                    continue

                if not same_dir:
                    weak_contigs = list(map(lambda c: c.reverse_copy(),
                                            weak_contigs[::-1]))
                link_to_change = left_scf.contigs[left_pos].link
                consistent.append(weak_contigs)

            if len(consistent) > 1:
                logger.debug("Something bad happening!")
            if len(consistent) != 1:
                continue

            new_contigs[-1].link = link_to_change
            new_contigs.extend(consistent[0])

        new_contigs.append(right_cnt)
        s = Scaffold(scf.name)
        s.contigs = new_contigs
        new_scafflods.append(s)

    return new_scafflods
