#!/usr/bin/env python

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


import source.graph_partiton as gp
import source.breakpoint_graph as bg
import source.sibelia_parser as sp
import source.overlap as ovlp
import source.debrujin_refine as debrujin
import source.scaffolder as scfldr
#from source.datatypes import Contig, Scaffold
from source.phylogeny import Phylogeny


def get_contigs(contigs_file):
    contigs = SeqIO.parse(contigs_file, "fasta")
    seqs = [seq for seq in contigs]
    names = [contig.id for contig in seqs]
    return seqs, names


def parse_config(filename):
    references = {}
    tree = None
    prefix = os.path.dirname(filename)

    for line in open(filename, "r"):
        line = line.strip()
        if not line:
            continue

        vals = line.split("=")
        if vals[0] == "tree":
            tree = vals[1].strip("\"")
        else:
            ref_id, filename = vals[0], vals[1].strip("\"")
            references[ref_id] = os.path.join(prefix, filename)

    return references, tree



def do_job(config_file, target_file, out_dir, block_size, skip_sibelia):
    if not os.path.isdir(out_dir):
        sys.stderr.write("Output directory doesn`t exists\n")
        return

    out_scaffolds = os.path.join(out_dir, "scaffolds.fasta")
    out_order = os.path.join(out_dir, "order.txt")
    out_ref_scaffolds = os.path.join(out_dir, "scaffolds_refined.fasta")
    out_ref_order = os.path.join(out_dir, "order_refined.txt")
    out_graph = os.path.join(out_dir, "breakpoint_graph.dot")
    out_overlap = os.path.join(out_dir, "contigs_overlap.dot")

    #debug_dir = os.path.join(out_dir, "debug")
    #if not os.path.isdir(debug_dir):
    #    os.mkdir(debug_dir)
    debug_dir = None

    references, tree_string = parse_config(config_file)

    contigs_seqs, contig_names = get_contigs(target_file)
    sibelia_output = sp.SibeliaRunner(references, target_file, block_size,
                                        out_dir, contig_names, skip_sibelia)

    graph = bg.BreakpointGraph()
    graph.build_from(sibelia_output)

    phylogeny = Phylogeny(tree_string)

    adj_finder = gp.AdjacencyFinder(graph, phylogeny, debug_dir)
    connections = adj_finder.find_adjacencies()
    scaffolds = scfldr.get_scaffolds(connections, sibelia_output)

    contigs_dict = {seq.id : seq.seq for seq in contigs_seqs}
    scfldr.output_scaffolds(contigs_dict, scaffolds, out_scaffolds, out_order)

    #ovlp.build_graph(contigs_dict, KMER, open(out_overlap, "w"))
    #refined_scaffolds = debrujin.refine_contigs(out_overlap, scaffolds)
    #output_scaffolds(contigs_dict, refined_scaffolds, out_ref_scaffolds,
    #                                            out_ref_order)


def main():
    parser = argparse.ArgumentParser(description="Tool for reference-assisted assembly")
    parser.add_argument("-c", action="store", metavar="config", dest="config",
                        required=True, help="Configuration file")
    parser.add_argument("-t", action="store", metavar="target_assembly", dest="target_assembly",
                        required=True, help="File with assembly in fasta format")
    parser.add_argument("-o", action="store", metavar="output_dir", dest="output_dir",
                        required=True, help="Output directory")
    parser.add_argument("-s", action="store_const", metavar="skip_sibelia", dest="skip_sibelia",
                        default=False, const=True, help="Skip Sibelia running step")
    parser.add_argument("-m", action="store", metavar="block_size", dest="block_size",
                        default="5000", help="Sibeia minimum block size (default: 5000)")
    #parser.add_argument("-g", action="store_const", metavar="debrujin_refine", dest="debrujin_refine",
    #                    default=False, const=True, help="Refine with Debrujin graph")

    args = parser.parse_args()
    do_job(args.config, args.target_assembly, args.output_dir, args.block_size, args.skip_sibelia)

if __name__ == "__main__":
    main()
