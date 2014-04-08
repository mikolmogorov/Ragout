Usage instructions for Ragout
=============================

    bin/ragout [-h] [-o OUTPUT_DIR] [-r] [-v] config_file
    
Supported arguments:

    positional arguments:
      config_file           path to the configuration file

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT_DIR, --outdir OUTPUT_DIR
                            path to the working directory (default: ragout-out)
      -s {sibelia,cactus}, --synteny {sibelia,cactus}
                            which tool to use for synteny block decomposition.
                            (default: sibelia)
      --refine              refine with the assembly graph (default: False)
      --circular            treat input references as circular (default: False)
      --overwrite           overwrite existing Sibelia/Cactus results (default:
                            False)
      --debug               enable debug output (default: False)
      --version             show program's version number and exit


Examples
---------

You can try Ragout on the provided ready-to-use examples:

    bin/ragout examples/E.Coli/ecoli.cfg --outdir examples/E.Coli/out/ --refine
    bin/ragout examples/H.Pylori/helicobacter.cfg --outdir examples/H.Pylori/out/ --refine
    bin/ragout examples/S.Aureus/aureus.cfg --outdir examples/S.Aureus/out/ --refine
    bin/ragout examples/V.Cholerea/cholerea.cfg --outdir examples/V.Cholerae/out/ --refine

Algorithm overview
------------------

This is a very brief description of the algorithm. See our paper 
for the detailed explanation.

Ragout works with genomes represented as sequences of synteny blocks
and firstly uses *Sibelia* for this decompostion. 
Next, Ragout assembles contigs into scaffolds using a breakpoint graph.

This procedure is repeated multiple times with the different size
of synteny block decomposition. Afterwards, an optional refinement
step is performed (if --refine was specified).

Input
------

Ragout takes as input:

- Reference sequences in *fasta* format
- Target assembly in *fasta* format (a set of contigs)
- Phylogenetic tree for both reference and target genomes in "newick" format
- Minimum synteny block size (in multiple scales)

All these parameters should be described in a single configuration file.
See the example of such file below.

Configuration file
------------------

Here is an example of Ragout configuration file:

    REF col=references/COL.fasta
    REF jkd=references/JKD6008.fasta
    REF rf122=references/RF122.fasta
    REF n315=references/N315.fasta

    TARGET usa=usa300_contigs.fasta

    TREE=(rf122:0.0280919,(((usa:0.0151257,col:0.0127906):0.0132464,jkd:0.0439743):0.00532819,n315:0.0150894):0.0150894);

    BLOCK=5000,500,100

Keywords description:

- REF: label of the reference sequence and its relative path
- TARGET: label of the target assembly and its relative path
- TREE: phylogenetic tree for both reference and target genomes
- BLOCK: minimum synteny block size (in multiple scales, one per iteration)


Output files
------------

After running Ragout, an output directory will contain:

* "scaffolds.ord" with a resulting order of contigs
* "scaffolds.fasta" with scaffold sequences (contigs are separated by 11 Ns)
* "scaffolds_refined.ord" with a refined order contigs (if --refine was specified)
* "scaffolds_refined.fasta" with refined scaffold sequences (if --refine was specified)

The parameters choice
---------------------

### Minimum synteny block size

Because the decomposition procedure is parameter-dependent, the assembly
is performed in multiple iterations with different synteny block
scale. Intuitively, the algorithm firstly considers only contigs
that are long enough and then puts shorter ones into the analysis.

For bacterial genomes, we recommend to run Ragout in three
iterations with the block size equal to 5000, 500, 100.
However, you can specify our own configuration which better
describes your dataset.

### Phylogenetic tree

Running with multiple references, the output of Ragout may highly
depend of the given phylogenetic tree and can be biased if
the tree is incorrect. We recommend to carefully assess the
tree before performing any significant studies.

In our future work we are planning to implement some assessment/correction
procedures for the phylogenetic tree.

Recovering synteny blocks from maf file
---------------------------------------

Ragout has a module which recovers synteny blocks from *MAF* multiple 
alignment file inside Ragout's pipeline. This module can also be run
as a standalone tool which allows you convert an arbitrary *MAF* file
into synteny blocks. Output format is similar to one from *Sibelia*.
To run this module, use:

    bin/maf2synteny maf_file output_dir synteny_block_size

Useful scripts
--------------

Scripts are located in "scripts" directory

**verify-order.py:**

Tests the correctness of the infered order of contigs if a closely related reference
is available. First, contigs should be mapped on this reference using *nucmer* software:

    nucmer --maxmatch --coords reference contigs

Then run the script with the obtained "coords" file:

    scripts/verify-order.py nucmer_coords ord_file