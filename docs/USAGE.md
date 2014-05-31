Usage instructions for Ragout
=============================

    Usage: ragout [-h] [-o OUTPUT_DIR] [-s {sibelia,cactus,maf}] [--refine]
                  [--overwrite] [--debug] [--version]
                  recipe_file
    
Supported arguments:

    positional arguments:
      recipe_file           path to recipe file

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT_DIR, --outdir OUTPUT_DIR
                            path to the working directory (default: ragout-out)
      -s {sibelia,cactus,maf}, --synteny {sibelia,cactus,maf}
                            which tool to use for synteny block decomposition.
                            (default: sibelia)
      --refine              refine with the assembly graph (default: False)
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

    .tree = (rf122:0.02,(((usa:0.01,col:0.01):0.01,jkd:0.04):0.005,n315:0.01):0.01);
    .target = usa
    .blocks = 5000, 500, 100

    *.circular = true

    col.fasta = references/COL.fasta
    jkd.fasta = references/JKD6008.fasta
    rf122.fasta = references/RF122.fasta
    n315.fasta = references/N315.fasta
    usa.fasta = usa300_contigs.fasta
    

or, if using *MAF* as input:

    .tree = (miranda:0.04,(yakuba:0.12,(melanogaster:0.04,simulans:0.04):0.08):0.15);
    .target = yakuba
    .maf = genomes/alignment.maf
    .blocks = 5000, 500

    yakuba.fasta = genomes/D.yakuba_contigs.fasta
   

###Parameters description:

Each parameter could be "global" or "local" (for a particular genome).
Global parameters starts from dot:

    .global_param_name = value

To set local parameter, use:

    genome_name.param_name = value

###Global parameters

* tree: phylogenetic tree in NEWICK format [required]
* target: target genome name [required]
* blocks: comma-separated list of minimum synteny block sizes [required]
* maf: path to multiple alignment in *MAF* format [default = not set]

###Local parameters

* fasta: path to *FASTA* with genomic sequences [default = not set]
* circular: indicates that reference chromosomes are circular [default = false]
* draft: indicated that reference is in draft form (contigs/scaffolds) [default = false]

###Default values

You can change default values for local parameters assigning parameter value to '*' genome.
For instance, if all input references except one are in draft form, you can write:

    *.draft = true
    complete_ref.draft = false

###Detailed description

Genome names are picked form the terminal nodes of the phylogenetic tree.
All those names should be uniqe. If the branch length is ommited, it would be set to 1.

Paths to *FASTA*/*MAF* can be both relative and absolute.

Ragout firstly decomposes genomes into set of synteny blocks.
You can use either a set of *FASTA* files corresponding to each input genome
or multiple alignment of all the genomes in *MAF* format.
In both cases you should specify target's *FASTA* since it will be
used to generate output. See "Synteny backends" section for more information.

Output files
------------

After running Ragout, an output directory will contain:

* "scaffolds.ord" with a resulting order of contigs
* "scaffolds.fasta" with scaffold sequences (contigs are separated by 11 Ns)
* "scaffolds_refined.ord" with a contigs order after refinement (if --refine was specified)
* "scaffolds_refined.fasta" with refined scaffold sequences (if --refine was specified)


Synteny backends
----------------

Ragout have three different options for synteny block decomposition:

* Decompoition with *Sibelia*
* Decomposition with *progressiveCactus*
* Use of finished alignment in *MAF* format

You can choose between backends by specifying --synteny (-s) option.
If you use *Sibelia* or *progressiveCactus*, you should specify separate *FASTA*
file for each input genome, while if you work with *MAF*, you should set
only a path to *MAF* itself and a path to targset's *FASTA* (see below).

### Sibelia

"Sibelia" option is set by default and is recommended for small genomes (like bacterial ones).
The tools should be installed in your system -- see docs/INSTALL.md for detailed instructions.

### progressiveCactus

"progressiveCactus" can be used for bigger genomes, up to multuple mammalian species.
Please note, that current implementation is still experimental. The tool also 
should be properly installed. Do not forget to mask repeats (with RepeatMasker, for instance)
before applying progressiveCactus to genomes with big fraction of repetitive sequences.

### alignment in *MAF* format

If you already have a multiple alignment, you also can use it for synteny blocks decomposition.
Alignment should be in *MAF* format and sequence names should follow UCSC notation:

    genome_name.sequence_name

In case you are working with *MAF* input you should not specify reference *FASTA* files.
All you need is to set *FASTA* for target genome (which will be used for output generation).


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

    .tree = (target,ref1:0.1,ref2:0.05,ref3:0.003);


Useful scripts
--------------

Scripts are located in "scripts" directory

**verify-order.py:**

Tests the correctness of the infered order of contigs if a closely related reference
is available. First, contigs should be mapped on this reference using *nucmer* software:

    nucmer --maxmatch --coords reference contigs

Then run the script with the obtained "coords" file:

    scripts/verify-order.py nucmer_coords ord_file