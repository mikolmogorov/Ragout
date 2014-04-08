#!/usr/bin/env python2.7

from __future__ import print_function
import sys, os
from collections import namedtuple, defaultdict
import networkx as nx
import pylab
from Bio import Phylo
from cStringIO import StringIO

from utils.nucmer_parser import *


def verify_alignment(alignment, contigs):
    problematic_contigs = []
    by_name = defaultdict(list)
    for entry in alignment:
        by_name[entry.contig_id].append(entry)
    for name in contigs:
        if len(by_name[name]) > 1:
            hits = list(map(lambda e: (e.s_ref, e.len_qry), by_name[name]))
            print("WARNING: Duplicated contig", name, hits, file=sys.stderr)
            problematic_contigs.append(name)
        if not by_name[name]:
            print("WARNING: Contig", name, "is not aligned", file=sys.stderr)
            problematic_contigs.append(name)
    return problematic_contigs


def get_true_adjacencies(alignment, contig_permutations, break_contigs, circular):
    by_chr = group_by_chr(alignment)
    adjacencies = []

    for chr_name, entries in by_chr.items():
        if circular:
            entries.append(entries[0])

        prev_block = None
        prev_contig = None
        for hit in entries:
            if prev_contig in break_contigs or hit.contig_id in break_contigs:
                continue

            sign = 1 if hit.e_qry > hit.s_qry else -1
            blocks = contig_permutations[hit.contig_id]
            #print(hit.contig_id, blocks)

            if sign < 0:
                blocks = list(map(lambda x: -x, blocks))[::-1]
            if prev_block:
                adjacencies.append((-prev_block, blocks[0]))
            prev_block = blocks[-1]
            prev_contig = hit.contig_id

    return adjacencies


def get_contig_permutations(filename):
    contigs = {}
    for line in open(filename, "r"):
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            name = line[1:]
        else:
            blocks = line.split(" ")[:-1]
            contigs[name] = list(map(int, blocks))
    return contigs


def output_edges(edges, out_file):
    fout = open(out_file, "w")
    fout.write("graph {\n")
    for (v1, v2) in edges:
        fout.write("{0} -- {1};\n".format(v1, v2))
    fout.write("}")


def g2c(genome_id):
    if genome_id not in g2c.table:
        g2c.table[genome_id] = g2c.colors[0]
        g2c.colors = g2c.colors[1:] + g2c.colors[:1] #rotate list
    return g2c.table[genome_id]
g2c.colors = ["green", "blue", "yellow", "cyan", "magnetta"]
g2c.table = {}


def draw_breakpoint_graph(base_dot, predicted_dot, true_edges, out_dir):
    base_graph = nx.read_dot(base_dot)
    predicted_edges = nx.read_dot(predicted_dot)
    out_graph = nx.MultiGraph()
    for v1, v2, data in base_graph.edges_iter(data=True):
        color = g2c(data["genome_id"])
        out_graph.add_edge(v1, v2, color=color)
    for v1, v2 in predicted_edges.edges_iter():
        out_graph.add_edge(v1, v2, color="red", style="dashed")
    for (v1, v2) in true_edges:
        out_graph.add_edge(str(v1), str(v2), color="red", style="bold")

    subgraphs = nx.connected_component_subgraphs(out_graph)
    for comp_id, subgr in enumerate(subgraphs):
        if len(subgr) == 2:
            continue
        comp_file = os.path.join(out_dir, "comp{0}-bg.png".format(comp_id))
        agraph = nx.to_agraph(subgr)
        agraph.layout(prog="dot")
        agraph.draw(comp_file)


def draw_phylogeny(phylogeny_txt, out_file):
    tree_string, target_name = open(phylogeny_txt, "r").read().splitlines()
    g2c.table[target_name] = "red"

    tree = Phylo.read(StringIO(tree_string), "newick")
    tree.clade.branch_length = 0
    for clade in tree.find_clades():
        if clade.is_terminal():
            clade.color = g2c(clade.name)
    tree.ladderize()
    pylab.rcParams["lines.linewidth"] = 3.0
    Phylo.draw(tree, do_show=False)

    pylab.savefig(out_file)


def do_job(nucmer_coords, debug_dir, circular):
    used_contigs = os.path.join(debug_dir, "used_contigs.txt")
    true_adj_out = os.path.join(debug_dir, "true_edges.dot")
    base_dot = os.path.join(debug_dir, "breakpoint_graph.dot")
    predicted_dot = os.path.join(debug_dir, "predicted_edges.dot")
    phylogeny_in = os.path.join(debug_dir, "phylogeny.txt")
    phylogeny_out = os.path.join(debug_dir, "phylogeny.png")

    draw_phylogeny(phylogeny_in, phylogeny_out)

    contigs = get_contig_permutations(used_contigs)
    alignment = parse_nucmer_coords(nucmer_coords)
    alignment = list(filter(lambda e: e.contig_id in contigs, alignment))
    #alignment = join_collinear(alignment)
    alignment = filter_by_coverage(alignment)
    alignment = join_collinear(alignment)
    break_contigs = verify_alignment(alignment, contigs)

    true_adj = get_true_adjacencies(alignment, contigs, break_contigs, circular)
    output_edges(true_adj, true_adj_out)
    draw_breakpoint_graph(base_dot, predicted_dot, true_adj, debug_dir)
    #print(g2c.table)


def main():
    if len(sys.argv) < 3:
        print("Usage: debug_report.py <nucmer_coords> <debug_dir>")
        return

    nucmer_coords = sys.argv[1]
    debug_dir = sys.argv[2]
    circular = True
    do_job(nucmer_coords, debug_dir, circular)

if __name__ == "__main__":
    main()
