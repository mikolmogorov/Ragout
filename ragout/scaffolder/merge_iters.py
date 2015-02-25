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

                if (left_scf != right_scf or
                        abs(left_pos - right_pos) == 1 or
                        big_sign != small_sign):
                    consistent = False

                same_dir = left_pos < right_pos
                if not same_dir:
                    left_pos, right_pos = right_pos, left_pos

                weak_contigs = left_scf.contigs[left_pos + 1 : right_pos]
                if any(c.seq_name in big_unique for c in weak_contigs):
                    consistent = False

                if not same_dir:
                    weak_contigs = list(map(lambda c: c.reverse_copy(),
                                            weak_contigs[::-1]))
                link_to_change = left_scf.contigs[left_pos].link

            new_contigs.append(left_cnt)
            if consistent:
                new_contigs[-1].link = link_to_change
                new_contigs.extend(weak_contigs)
                #logger.debug("Inserting '{0}' between {1} and {2}"
                #             .format(map(lambda c: c.seq_name, weak_contigs),
                #                     left_cnt, right_cnt))
            else:
                new_contigs.extend(big_scf.contigs[left_idx+1:right_idx])

        new_contigs.append(right_cnt)
        s = Scaffold(big_scf.name)
        s.contigs = new_contigs
        new_scafflods.append(s)

    return new_scafflods
