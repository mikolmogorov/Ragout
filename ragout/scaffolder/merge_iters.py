#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module merges scaffolds from two consecutive
iterations
"""

from collections import namedtuple, defaultdict
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
            big_index.add(c.name)

    small_index = {}
    for scf in small_scaffolds:
        for pos, contig in enumerate(scf.contigs):
            assert contig.name not in small_index
            small_index[contig.name] = (scf, pos)

    count = 0
    new_scafflods = []
    for scf in big_scaffolds:
        result = []
        for prev_cont, new_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            result.append(prev_cont)

            try:
                scf_prev, begin = small_index[prev_cont.name]
                scf_new, end = small_index[new_cont.name]
            except KeyError:
                continue
            if scf_prev.name != scf_new.name:
                continue

            assert end != begin
            same_dir = True
            if end < begin:
                same_dir = False
                end, begin = begin, end

            consistent = True
            for c in scf_prev.contigs[begin + 1 : end]:
                if c.name in big_index:
                    consistent = False
                    break

            if not consistent or end - begin == 1:
                continue

            if ((prev_cont.sign == new_cont.sign) !=
                (scf_prev.contigs[begin].sign == scf_prev.contigs[end].sign)):
                continue

            count += end - begin - 1
            contigs = scf_prev.contigs[begin + 1 : end]
            if not same_dir:
                contigs = contigs[::-1]
                contigs = list(map(lambda c: c.reverse(), contigs))
            #keeping gap from new contigs
            result[-1].link = scf_prev.contigs[begin].link
            result.extend(contigs)

        result.append(new_cont)
        s = Scaffold(scf.name)
        s.contigs = result
        new_scafflods.append(s)

    return new_scafflods
