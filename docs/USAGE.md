Usage instructions for Ragout
=============================

    Usage: ragout [-h] [-o OUTPUT_DIR] [-s {sibelia,cactus}] [--refine]
                  [--circular] [--overwrite] [--debug] [--version]
                  recipe_file
    
Supported arguments:

    positional arguments:
      recipe_file           path to recipe file

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

    ragout examples/E.Coli/ecoli.rcp --outdir examples/E.Coli/out/ --refine
    ragout examples/H.Pylori/helicobacter.rcp --outdir examples/H.Pylori/out/ --refine
    ragout examples/S.Aureus/aureus.rcp --outdir examples/S.Aureus/out/ --refine
    ragout examples/V.Cholerea/cholerea.rcp --outdir examples/V.Cholerae/out/ --refine

Algorithm overview
------------------

This is a very brief description of the algorithm. See our paper 
for the detailed explanation.

Ragout works with genomes represented as sequences of synteny blocks
and firstly uses *Sibelia* for this decomposition. 
Next, Ragout assembles contigs into scaffolds using a breakpoint graph.

This procedure is repeated multiple times with the different size
of synteny block decomposition. Afterwards, an optional refinement
step with assembly (overlap) graph is performed (if --refine was specified).

Input
------

Ragout takes as input:

* Reference genomes in *FASTA* format
* Target (assembling) genome in *FASTA* format (a set of contigs)
* Phylogenetic tree for both reference and target genomes in *NEWICK* format
* Minimum synteny block size (in multiple scales)

All these parameters should be described in a single recipe file.
See the example of such file below.

Recipe file
-----------

If you want to cook Ragout, you need to write a recipe first.
Here is an example of such recipe file:

    REF col = references/COL.fasta
    REF jkd = references/JKD6008.fasta
    REF rf122 = references/RF122.fasta
    REF n315 = references/N315.fasta

    TARGET usa = usa300_contigs.fasta

    TREE = (rf122:0.0280919,(((usa:0.0151257,col:0.0127906):0.0132464,jkd:0.0439743):0.00532819,n315:0.0150894):0.0150894);

    BLOCK = 5000,500,100

Keywords description:

* REF: reference genome name and a path to FASTA file with sequences
* TARGET: target genome name and a path to FASTA file with sequences
* TREE: phylogenetic tree in NEWICK format
* BLOCK: minimum synteny block sizes (in multiple scales, one per iteration)

All names should be uniqe. Paths can be both relative and absolute.
The tree should contain all the described genomes (both references and target)
in leaf nodes. If the branch length is ommited, it would be set to 1.
NEWICK string as well as BLOCK string should not contain any spaces inside.

Output files
------------

After running Ragout, an output directory will contain:

* "scaffolds.ord" with a resulting order of contigs
* "scaffolds.fasta" with scaffold sequences (contigs are separated by 11 Ns)
* "scaffolds_refined.ord" with a contigs order after refinement (if --refine was specified)
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
the tree is incorrect.

If the phylogeny is unknown or ambiguous, you are still able run Ragout assuming
the "star" phylogeny and specifying the evolutionary distance between
target and references (which is easier to find out):

    TREE = (target,ref1:0.1,ref2:0.05,ref3:0.003);


Experimental support of Progressive Cactus
------------------------------------------
As Sibelia was designed for bacterial species comparison, we are planning
to use Progressive Cactus for recovering synteny blocks in bigger genomes.
You already can try it, however this option is still quite unstable.
Note, that repeats should be masked before running Ragout on big genomes
(with RepeatMasker, for instance).

First, download and build Progressive Cactus: https://github.com/glennhickey/progressiveCactus
Then run Ragout with CACTUS_INSTALL pointing to cactus installation directory:

    export CACTUS_INSTALL=your_cactus_dir
    ragout -s cactus ...


Useful scripts
--------------

Scripts are located in "scripts" directory

**verify-order.py:**

Tests the correctness of the infered order of contigs if a closely related reference
is available. First, contigs should be mapped on this reference using *nucmer* software:

    nucmer --maxmatch --coords reference contigs

Then run the script with the obtained "coords" file:

    scripts/verify-order.py nucmer_coords ord_file