Usage instructions for Ragout
=============================

    python ragout.py [-h] -c config -o output_dir [-s] [-g]
    
Supported arguments:

    -h, --help      show the help message and exit
    -c config       Configuration file
    -o output_dir   Output directory
    -s              Skip Sibelia running step
    -g              Refine with assembly graph

Examples
---------

You can try Ragout on the provided ready-to-use examples:

    python ragout.py -c examples/E.Coli/ecoli.cfg -o examples/E.Coli/out/
    python ragout.py -c examples/H.Pylori/helicobacter.cfg -o examples/H.Pylori/out/
    python ragout.py -c examples/S.Aureus/aureus.cfg -o examples/S.Aureus/out/
    python ragout.py -c examples/V.Cholerea/cholerea.cfg -o examples/V.Cholerea/out/

Algorithm overview
------------------

This is a very brief description of the algorithm. See our paper 
for the detailed explanation.

Ragout works with genomes represented as sequences of synteny blocks
and firstly uses *Sibelia* for this decompostion. 
Next, Ragout assembles contigs into scaffolds using a breakpoint graph.

This procedure is repeated multiple times with the different size
of synteny block decomposition. Afterwards, an optional refinement
step is performed (is -g is specified).

Input
------

Ragout takes as input:

- Reference sequences in *fasta* format
- Target assembly in *fasta* format (a set of contigs)
- Phylogenetic tree for both reference and target genomes in "newick" format
- Minimum synteny block size (in multiple scales)

All these settings should be described in a single config file.
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

Keywords explanation:

- REF: label of the reference sequence and its relative path
- TARGET: label of the target assembly and its relative path
- TREE: phylogenetic tree for both reference and target genomes
- BLOCK: minimum synteny block size (in multiple scales, one per iteration)


Output files
------------

After running Ragout, an output directory will contain:

* "scaffolds.ord" with a resulting order of contigs
* "scaffolds.fasta" with scaffold sequences (contigs are separated by 11 Ns)
* "scaffolds_refined.ord" with a refined order contigs (if you used -g option)
* "scaffolds_refined.fasta" with refined scaffold sequences (if you used -g option)

The parameters choice
---------------------

### Minimum synteny block size

Because the decomposition procedure is parameter-dependent 
(it requires a certain minimum synteny block size), the assembly
is performed in multiple iterations with different synteny block
scale. Intuitively, the algorithm firstly considers only contigs
that are long enough and then puts shorter ones into the analysis.

For bacterial genomes, we recommend to run Ragout in three
iterations with the block size equal to 5000, 500, 100.
However, you can specify our own configuration which better
describes your dataset.

### Phylogenetic tree

With multiple references, the algorithm's output may be highly 
dependent of the phylogenetic tree structure and an incorrect
choice of it can bias the result. We recommend to carefully assess the
tree before doing any significant studies.

In our future work we are planning to implement some assessment/correction
procedures for the phylogenetic tree.

Useful scripts
--------------

Scripts are located in "scripts" directory

**test.py:**

Tests the correctness of the infered contigs order if a closely related reference
is available. First, contigs should be mapped on this reference using *nucmer*:

    nucmer --maxmatch --coords reference contigs

Then run the script with the obtained "coords" file:

	python test.py nucmer_coords ord_file