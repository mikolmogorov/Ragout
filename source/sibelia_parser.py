import os
import shutil
import subprocess
import copy
from collections import namedtuple, defaultdict
from Bio import SeqIO


SIBELIA_BIN = "/home/volrath/Bioinf/Sibelia/distr/bin/Sibelia"


def run_sibelia(fasta_files, block_size, out_dir):
    print "Running Sibelia..."
    cmdline = [SIBELIA_BIN, "-s", "loose", "-m", str(block_size), "-o", out_dir]
    cmdline.extend(fasta_files)
    process = subprocess.Popen(cmdline)
    process.wait()
    os.remove(os.path.join(out_dir, "coverage_report.txt"))
    os.remove(os.path.join(out_dir, "d3_blocks_diagram.html"))
    os.remove(os.path.join(out_dir, "blocks_coords.txt"))
    shutil.rmtree(os.path.join(out_dir, "circos"))
    return os.path.join(out_dir, "genomes_permutations.txt")


def get_chr_names(genomes):
    chr_to_id = {}
    for seq_id, seq_file in genomes.iteritems():
        for seq in SeqIO.parse(seq_file, "fasta"):
            chr_to_id[seq.id] = seq_id
    return chr_to_id


def split_permutations(chr_to_gen, references, targets, perm_file, out_dir):
    out_files = {}
    config = open(os.path.join(out_dir, "blocks.cfg"), "w")
    genomes = dict(references.items() + targets.items())

    for gen_id in set(chr_to_gen.values()):
        filename = genomes[gen_id]
        base = os.path.splitext(os.path.basename(filename))[0]
        block_file_base = base + ".blocks"
        block_file = os.path.join(out_dir, block_file_base)

        out_files[gen_id] = open(block_file, "w")
        if gen_id in references:
            config.write("REF {0}={1}\n".format(gen_id, block_file_base))
        else:
            assert gen_id in targets
            config.write("TARGET {0}={1}\n".format(gen_id, block_file_base))

    for line in open(perm_file, "r").read().splitlines():
        if line.startswith(">"):
            name = line[1:]
        else:
            handle = out_files[chr_to_gen[name]]
            handle.write(">{0}\n{1}\n".format(name, line))



def make_permutations(references, targets, block_size, output_dir):
    genomes = dict(references.items() + targets.items())
    perm_file = run_sibelia(genomes.values(), block_size, output_dir)
    #perm_file = os.path.join(output_dir, "genomes_permutations.txt")
    chr_to_gen = get_chr_names(genomes)
    split_permutations(chr_to_gen, references, targets, perm_file, output_dir)

