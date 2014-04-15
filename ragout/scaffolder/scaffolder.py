#This module assembles contigs into scaffolds with respect
#to given adjacencies. Also, it outputs scaffolds in different
#formats
#################################################################

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import copy
import logging

from ragout.shared.datatypes import Contig, Scaffold

logger = logging.getLogger()

#PUBLIC:
################################################################

#Assembles scaffolds
def get_scaffolds(connections, perm_container):
    logger.info("Building scaffolds")
    scaffolds = _extend_scaffolds(connections, perm_container)
    scaffolds = list(filter(lambda s: len(s.contigs) > 1, scaffolds))
    return scaffolds


#Outputs scaffolds to file in "ord" format
def output_order(scaffolds, out_order):
    out_order_stream = open(out_order, "w")
    for scf in scaffolds:
        out_order_stream.write(">" + scf.name + "\n")
        for contig in scf.contigs:
            out_order_stream.write(str(contig) + "\n")


#Outputs scaffodls to file in "fasta" format
def output_fasta(target_dict, scaffolds, out_fasta):
    MIN_CONTIG_LEN = 0

    contigs_fasta = {}
    for target_file in target_dict.values():
        for seq in SeqIO.parse(target_file, "fasta"):
            contigs_fasta[seq.id] = seq.seq

    out_stream = open(out_fasta, "w")
    used_contigs = set()

    for scf in scaffolds:
        scf_seq = Seq("")
        first = True

        for contig in scf.contigs:
            cont_seq = contigs_fasta[contig.name]
            used_contigs.add(contig.name)

            if contig.sign < 0:
                cont_seq = cont_seq.reverse_complement()

            if not first:
                scf_seq += Seq("N" * 11)
            first = False
            scf_seq += cont_seq

        SeqIO.write(SeqRecord(scf_seq, id=scf.name, description=""), out_stream, "fasta")

    count = 0
    for h, seq in contigs_fasta.items():
        if len(seq) > MIN_CONTIG_LEN and h not in used_contigs:
            count += 1


#PRIVATE:
################################################################

#Assembles contigs into scaffolds
def _extend_scaffolds(connections, perm_container):
    contigs, contig_index = _make_contigs(perm_container)

    scaffolds = []
    visited = set()
    counter = [0]

    def extend_scaffold(contig):
        visited.add(contig)
        scf_name = "scaffold{0}".format(counter[0])
        counter[0] += 1
        scf = Scaffold.with_contigs(scf_name, contig.blocks[0], contig.blocks[-1], [contig])
        scaffolds.append(scf)

        #go right
        while scf.right in connections:
            adjacent = connections[scf.right].end

            #print adjacent, contig_index[abs(adjacent)]
            assert len(contig_index[abs(adjacent)]) == 1

            contig = contig_index[abs(adjacent)][0]
            if contig in visited:
                break

            if contig.blocks[0] == adjacent:
                scf.contigs.append(contig)
                scf.right = contig.blocks[-1]
                visited.add(contig)
                continue

            if -contig.blocks[-1] == adjacent:
                scf.contigs.append(contig)
                scf.contigs[-1].sign = -1
                scf.right = -contig.blocks[0]
                visited.add(contig)
                continue

            break

        #go left
        while -scf.left in connections:
            adjacent = -connections[-scf.left].end

            assert len(contig_index[abs(adjacent)]) == 1

            contig = contig_index[abs(adjacent)][0]
            if contig in visited:
                break

            if contig.blocks[-1] == adjacent:
                scf.contigs.insert(0, contig)
                scf.left = contig.blocks[0]
                visited.add(contig)
                continue

            if -contig.blocks[0] == adjacent:
                scf.contigs.insert(0, contig)
                scf.contigs[0].sign = -1
                scf.left = -contig.blocks[-1]
                visited.add(contig)
                continue

            break

    for contig in contigs:
        if contig not in visited:
            extend_scaffold(contig)
    return scaffolds


#Converts permutations into contigs
def _make_contigs(perm_container):
    contigs = []
    index = defaultdict(list)
    for perm in perm_container.target_perms_filtered:
        if len(perm.blocks) == 0:
            continue

        contigs.append(Contig(perm.chr_id))
        contigs[-1].blocks = copy.copy(perm.blocks)

        for block in perm.blocks:
            index[abs(block)].append(contigs[-1])

    return contigs, index
