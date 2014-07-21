Usage instructions for Ragout
=============================

    Usage: ragout.py [-h] [-o OUTPUT_DIR] [-s {sibelia,cactus,maf}] [--refine]
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

    python ragout.py examples/E.Coli/ecoli.rcp --outdir examples/E.Coli/out/ --refine
    python ragout.py examples/H.Pylori/helicobacter.rcp --outdir examples/H.Pylori/out/ --refine
    python ragout.py examples/S.Aureus/aureus.rcp --outdir examples/S.Aureus/out/ --refine
    python ragout.py examples/V.Cholerea/cholerea.rcp --outdir examples/V.Cholerae/out/ --refine


Sequence data
-------------

Ragout takes as input contigs from a short-read assembly. We performed
our tests with SPAdes, ABySS and SOAPdenovo assemblers, others should
work fine too, if their output satisfy the following conditions:

* Assembly coverage should be acceptable (80-90%+) and contigs/scaffolds
  should not overlap (except by ends, see below)
* *All* contigs/scaffolds output by assembler should be used (including
  short ones)
* For the better performance of the refinment module, contigs/scaffolds
  that were connected in a graph used by assembler should overlap on a
  certain constant value (usully k-mer or (k-1)-mer). This holds true
  for the most of assemblers which utilize de Bruijn graphs. Currently,
  for other types of assemblers (such as SGA) support of
  refinement procedure is limited.


Algorithm overview
------------------

This is a very brief description of the algorithm. See our paper 
for the detailed explanation.

Ragout works with genomes represented as sequences of synteny blocks
and firstly uses *Sibelia* or local multiple sequence alignment
for this decomposition. 

Next, Ragout constructs a breakpoint graph and predicts missing 
adjacencies between synteny blocks in target genome (as it is fragmented
into contigs, some adjacencies are missing). Then, contigs are being 
extended into scaffolds with respect to inferred adjacencies. This 
procedure is repeated multiple times with the different scale of synteny
block decomposition. 

Afterwards, a refinement step is performed. Ragout utilize assembly (overlap) 
graph which is reconstructed by overlapping input contigs/scaffolds. This 
helps to insert very short/repetitive contigs into the assembly.


Input
-----

Ragout takes as input:

* Reference genomes in *FASTA* format
* Target genome in *FASTA* format (a set of contigs/scaffolds)
* Phylogenetic tree containing reference and target genomes in *NEWICK* format
* Synteny block sizes (in multiple scales)

All these parameters should be described in a single recipe file.
See the example of such file below.


Output
------

After running Ragout, an output directory will contain:

* __scaffolds.ord__: contigs order
* __scaffolds.fasta__: scaffolds sequences
* __scaffolds_refined.ord__: contigs order after refinement (if --refine was specified)
* __scaffolds_refined.fasta__: refined scaffolds sequences (if --refine was specified)


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
Global parameters start from dot:

    .global_param_name = value

To set local parameter, use:

    genome_name.param_name = value

###Global parameters

* __tree__: phylogenetic tree in NEWICK format [required]
* __target__: target genome name [required]
* __blocks__: comma-separated list of minimum synteny block sizes [required]
* __maf__: path to local multiple sequence alignment in *MAF* format [default = not set]

###Local parameters

* __fasta__: path to *FASTA* with genomic sequences [default = not set]
* __circular__: indicates that reference chromosomes are circular [default = false]
* __draft__: indicates that reference is in draft form (contigs/scaffolds) [default = false]

###Default values

You can change default values for local parameters by assigning the 
parameter value to the special "star" object.
For instance, if all input references except one are in draft form, you can write:

    *.draft = true
    complete_ref.draft = false

###Detailed description

Genomes are picked form the terminal nodes of the phylogenetic tree.
All those names should be uniqe. If the branch length is ommited, it would
be set to 1.

Paths to *FASTA*/*MAF* can be both relative and absolute. Running with 
Sibelia requires all sequence headers among ALL *FASTA* files to be unique.

Ragout firstly decomposes genomes into set of synteny blocks.
You can use either a set of *FASTA* files corresponding to each input genome
or local multiple alignment of all the genomes in *MAF* format.
In both cases you should specify target's *FASTA* since it will be
used to generate output. See "Synteny backends" section for more information.


The parameters choice
---------------------

### Minimum synteny block size

Because the decomposition procedure is parameter-dependent, the assembly
is performed in multiple iterations with different synteny block
scale. Intuitively, the algorithm firstly considers only contigs
that are long enough and then insert shorter ones in the assembly.

For bacterial genomes, we recommend to run Ragout in three
iterations with the block size equal to 5000, 500, 100.
However, you can specify our own configuration which better
describes your dataset.

### Phylogenetic tree

Running with multiple references, the output of Ragout may highly
depend of the given phylogenetic tree and can be biased if
the tree is incorrect.

If the phylogeny is unknown or ambiguous, you are still able run 
Ragout assuming the "star" phylogeny and specifying the evolutionary
distance between target and references (which is easier to find out):

    .tree = (target, ref1:0.1, ref2:0.05, ref3:0.003);

### Circular genomes

If you are working with circular genomes (like bacterial ones) it is 
recommended to set corresponding parameter in recipe file (see previous
seqction). This would generate some extra adjacencies which could 
be useful.

### Reference genome in draft form

Ragout can use even a draft assembly (contigs/scaffolds) as reference 
sequences. If so, you should specify it in recipe file (see previous 
section).


Synteny backends
----------------

Ragout have three different options for synteny block decomposition:

* Decompoition with *Sibelia*
* Decomposition with *progressiveCactus*
* Use of external local multiple sequence alignment (in *MAF* format)

You can choose between backends by specifying --synteny (-s) option.
If you use *Sibelia* or *progressiveCactus*, you should specify separate
*FASTA* file for each input genome, while if you work with *MAF*, you
should set only a path to *MAF* itself and a path to targset's *FASTA*
(see below).

### Sibelia

"Sibelia" option is set by default and is recommended for small genomes
(like bacterial ones).

### progressiveCactus

"progressiveCactus" can be used for bigger genomes, up to multuple 
mammalian species. Please note, that current implementation is still 
experimental. The tool also should be properly installed. Do not forget 
to mask repeats (with RepeatMasker, for instance) before applying 
*progressiveCactus* to genomes with a big fraction of repetitive sequences.

### Local alignment in *MAF* format

If you already have a local multiple alignment, you also can use it for
synteny blocks decomposition. Alignment should be in *MAF* format and 
sequence names should follow UCSC notation:

    genome_name.sequence_name

In case you are working with *MAF* input you should not specify *FASTA*
files for references. All you need is to set *FASTA* for target genome 
(which will be used for output generation and refinement).


Refinement with assembly graph
------------------------------

Ragout uses assembly (overlap) graph to incorporate very short / repetitive
contigs into assembly. First, this graph is reconstructed by overlapping
input contigs/scaffolds (see Sequence data seqction for requirements).
Then ragout scaffolds which are already available are being "threaded"
through this graph to find the true "genome path".

This procedure:

* Increases number of conitgs in output scaffolds
* Improves estimation of distances between contigs (will be filled with Ns)

However, sometimes the assembly graph is not accurate: some adjacencies
between contigs could be missing and on the other hand, there also might
be some false-positive adjacencies. This may lead to some incorrectly inserted
contigs. 

For good assemblies (with reasonable coverage, read length etc.) the fraction
of such errors should be very small or even zero. And even they exist,
they are "local" and do not violate the genome structure (probably, most
of them even will not be detected as missassembles by quality assesment
tools like Quast).


Useful scripts
--------------

Scripts are located in "scripts" directory

**verify-order.py:**

Tests the correctness of the infered contigs order if a "true" reference
is available. First, contigs should be mapped on that reference using 
*nucmer* software:

    nucmer --maxmatch --coords reference contigs

Then run the script with the obtained "coords" file:

    scripts/verify-order.py nucmer_coords ord_file