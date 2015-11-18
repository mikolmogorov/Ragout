#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)


from collections import defaultdict

def reestimate_gaps(scaffolds, bp_graph, phylogeny, fragment_seqs):
    inter_block = 0
    for scf in scaffolds:
        for ctg in scf.contigs:
            for b_1, b_2 in zip(ctg.perm.blocks[:-1], ctg.perm.blocks[1:]):
                v_1, v_2 = -b_1.signed_id(), b_2.signed_id()

                genomes = bp_graph.genomes_support(v_1, v_2)
                if len(genomes) > 1:
                    graph_dist = bp_graph.get_distance(v_1, v_2, phylogeny)
                    perm_dist = b_2.start - b_1.end

                    inter_block += graph_dist - perm_dist

    print(inter_block)


def check_blocks(permutation_container):
    blocks_diff = defaultdict(int)
    for perm in permutation_container.ref_perms:
        if perm.genome_name == "C57B6J":
            for bl in perm.blocks:
                blocks_diff[bl.block_id] += bl.length()

    for perm in permutation_container.target_perms:
        for bl in perm.blocks:
            blocks_diff[bl.block_id] -= bl.length()

    print(sum(blocks_diff.values()))
