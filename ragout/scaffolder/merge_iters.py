#This module merges scaffolds from two consecutive
#iterations
#########################################################

from collections import namedtuple
import sys
import logging

from ragout.shared.datatypes import Contig, Scaffold

logger = logging.getLogger()

#PUBLIC:
########################################################

#The only function here
def merge(big_scaffolds, small_scaffolds):
    logger.info("Merging two iterations")
    big_index = set()
    for scf in big_scaffolds:
        for c in scf.contigs:
            big_index.add(c.name)

    count = 0
    new_scafflods = []
    for scf in big_scaffolds:
        result = []
        for prev_cont, new_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            found_pair = False
            for small_scf in small_scaffolds:
                names = list(map(lambda c: c.name, small_scf.contigs))
                try:
                    begin = names.index(prev_cont.name)
                    end = names.index(new_cont.name)
                    found_pair = True
                    break
                except ValueError:
                    continue

            result.append(prev_cont)
            if not found_pair:
                continue

            assert end != begin
            same_dir = True
            if end < begin:
                same_dir = False
                end, begin = begin, end

            consistent = True
            for c in small_scf.contigs[begin + 1 : end]:
                if c.name in big_index:
                    consistent = False
                    break

            if not consistent or end - begin == 1:
                continue

            if ((prev_cont.sign == new_cont.sign) !=
                (small_scf.contigs[begin].sign == small_scf.contigs[end].sign)):
                continue

            count += end - begin - 1
            contigs = small_scf.contigs[begin + 1 : end]
            if not same_dir:
                contigs = contigs[::-1]
                contigs = list(map(lambda c: Contig(c.name, -c.sign, 0),
                                   contigs))
            result.extend(contigs)

        result.append(new_cont)
        s = Scaffold(scf.name)
        s.contigs = result
        new_scafflods.append(s)

    return new_scafflods
