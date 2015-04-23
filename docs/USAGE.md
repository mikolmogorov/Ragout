Usage Instructions for Ragout
=============================

Quick Usage
-----------

    usage: ragout.py [-h] [-o output_dir] [-s {sibelia,hal}]
                     [--no-refine] [--overwrite] [--repeats] [--debug]
                     [-t THREADS] [--version]
                     recipe_file


    positional arguments:
      recipe_file           path to recipe file

    optional arguments:
      -h, --help            show this help message and exit

      -o output_dir, --outdir output_dir
                            path to the working directory (default: ragout-out)
    
      -s {sibelia,hal}, --synteny {sibelia,hal}
                            backend for synteny block decomposition (default:
                            sibelia)
    
      --no-refine           disable refinement with assembly graph (default:
                            False)
    
      --overwrite           overwrite existing synteny blocks (default: False)
    
      --repeats             try to resolve repeats before constructing breakpoint
                            graph (default: False)
    
      --debug               enable debug output (default: False)
    
      -t THREADS, --threads THREADS
                            number of threads for synteny backend (default: 1)
    
      --version             show program's version number and exit


Examples
--------

You can try Ragout on the provided ready-to-use examples:

    python ragout.py examples/E.Coli/ecoli.rcp --outdir examples/E.Coli/out/
    python ragout.py examples/H.Pylori/helicobacter.rcp --outdir examples/H.Pylori/out/
    python ragout.py examples/S.Aureus/aureus.rcp --outdir examples/S.Aureus/out/
    python ragout.py examples/V.Cholerae/cholerae.rcp --outdir examples/V.Cholerae/out/


Assembly Requirements
---------------------

Ragout takes assembly fragments (contigs/scaffolds) as input. We performed
our tests with SPAdes, ABySS and SOAPdenovo assemblers, others should
work fine too, if their output satisfy the following conditions:

* Assembly coverage should be sufficient (80-90%+)
* *All* contigs/scaffolds output by assembler should be used (including
  short ones)
* For the better performance of the refinment module, contigs/scaffolds
  that were connected in a graph used by assembler should overlap on a
  constant value (usually k-mer or (k-1)-mer). This works
  for the most of assemblers which utilize de Bruijn graphs. Currently,
  for other types of assemblers (such as SGA) the performance of
  the refinement procedure is limited.

One current limitation of Ragout is that the algorithm is sensitive to
chimeric sequences. They may be expanded, which will result in 
artifically fused scaffolds. We recommend to reduce the fraction of chimeric
sequences as much as possible before running Ragout.



Algorithm Overview
-------------------

Ragout works with genomes represented as sequences of synteny blocks
and firstly uses *Sibelia* or *HAL alignment* for this decomposition. 
This step is usually the most time-consuming.

Next, Ragout constructs a breakpoint graph and predicts missing 
adjacencies between synteny blocks in target genome (as it is fragmented, 
some adjacencies are missing). Then, assembly fragments are being 
extended into scaffolds with respect to the inferred adjacencies. This 
procedure is repeated multiple times with the different scale of the synteny
blocks decomposition. 

Afterwards, a refinement step is performed. Ragout reconstruct 
assembly (overlap) graph by overlapping input assembly fragments. This 
graph is used to insert very short/repetitive fragments into the 
final scaffolds.


Input
-----

Ragout takes as input:

* Reference genomes [in *FASTA* format or packed into *HAL*]
* Target assembly in [in *FASTA* format or packed into *HAL*]

Optionally, you can add:

* Phylogenetic tree containing reference and target genomes [in *NEWICK* format]
* Synteny block scale

All these parameters should be described in a single recipe file.
See the example of such file below.


Output
------

After running Ragout, output directory will contain:

* __scaffolds.links__: scaffolds description (see below)
* __scaffolds.fasta__: scaffolds sequences
* __scaffolds_refined.links__: scaffolds description after refinement (see below)
* __scaffolds_refined.fasta__: scaffolds sequences after refinement



Recipe File
-----------

If you want to cook Ragout, you need to write a recipe first.
Here is an example of such recipe file:

    .references = rf122,col,jkd,n315
    .target = usa

    col.fasta = references/COL.fasta
    jkd.fasta = references/JKD6008.fasta
    rf122.fasta = references/RF122.fasta
    n315.fasta = references/N315.fasta
    usa.fasta = usa300_contigs.fasta

    *.circular = true
    .tree = (rf122:0.02,(((usa:0.01,col:0.01):0.01,jkd:0.04):0.005,n315:0.01):0.01);
    .blocks = small
    

or, if using *HAL* as input, tree and blocks scale are inffered automatically

    .references = miranda,simulans,melanogaster
    .target = yakuba
    .hal = genomes/alignment.hal
   

###Parameters description:

Each parameter could be "global" or "local" (for a particular genome).
Global parameters start from dot:

    .global_param_name = value

To set local parameter, use:

    genome_name.param_name = value

###Global parameters

* __references__: comma-separated list of reference names [*required*]
* __target__: target genome name [*required*]
* __tree__: phylogenetic tree in NEWICK format
* __blocks__: syntany block scale
* __hal__: path to the alignment in *HAL* format

###Local parameters

* __fasta__: path to *FASTA* [default = not set]
* __circular__: indicates that reference chromosomes are circular [default = false]
* __draft__: indicates that reference is in draft form (contigs/scaffolds) [default = false]

###Default values

You can change default values for local parameters by assigning the 
parameter value to the special "star" object.
For instance, if all input references except one are in draft form, you can write:

    *.draft = true
    complete_ref.draft = false

###Quick comments

Paths to *FASTA*/*HAL* can be both relative and absolute. 

If you use *Sibelia* or *Progressive Cactus*, you must specify 
FASTA for each input genome. If you use *HAL*, everything is taken from it.

Running with Sibelia requires all sequence headers (">gi...") 
among ALL *FASTA* files to be unique.

If you do not specify phylogenetic tree or synteny block scale, 
they are inferred automatically.


The Parameters Choice
---------------------

### Synteny block size

Because the decomposition procedure is parameter-dependent, the assembly
is performed in multiple iterations with different synteny block
scale. Intuitively, the algorithm firstly considers only fragments
that are long enough and then insert shorter ones into final scaffolds.

There are two pre-defined scales: "small" and "large". We recommend
"small" for relatively small genomes (bacterial) and large otherwise 
(mammalian). If th parameter is not set, it is automatically inferred
based on input genomes size.


### Phylogenetic tree

If this parameter is omitted, Ragout infers phylogenetic tree 
based on breakpoints data. Generally, it gives a good approximation.
However, if you already know the tree, we recommended to
guide the algorithm with it.

### Circular genomes

If you are working with circular genomes (like bacterial ones) it is 
recommended to set corresponding parameter in recipe file (see previous
section).

### Reference genome in draft form

Ragout can use even a draft assembly (contigs/scaffolds) as reference 
sequences. If so, you should specify it in recipe file (see previous 
section).


Synteny backends
----------------

Ragout have two different options for synteny block decomposition:

* Decomposition with *Sibelia*
* Use of whole genome alignment (in *HAL* format)

You can choose between backends by specifying --synteny (-s) option.
If you use *Sibelia*, you should specify separate *FASTA* file for 
each input genome, while if you work with *HAL*, it is not necessary.

### Sibelia

"Sibelia" option is set by default and is recommended for bacterial genomes.

### Whole genome alignment in *HAL* format

Alternatively, Ragout can use *HAL* whole genome alignment for synteny blocks
decomposition. This is recommended for large (mammalian) genomes, which
*Sibelia* can not process. This alignment is ouput by Progressive Cactus 
aligner [https://github.com/glennhickey/progressiveCactus].

The pipeline is as follows: first, align references and target genomes using
Progressive Cactus (you will need phylogenetic tree) and then run Ragout
on the resulted *HAL* file. This would require *HAL Tools* package to be 
installed in your system.

### Progressive Cactus and MAF backends are deprecated

The support of MAF synteny backend is deprecated, because it is more 
convenient to work directly with *HAL*, which is a default output of
Progressive Cactus.

Progressive Cactus backend (not to be confused with HAL backend) 
is also deprecated, because typical run of the aligner
includes some parallelization steps, which are hard to automate.

Please let us know, if you have any issues with that changes.

Repeat Resolution
-----------------
As the main Ragout algorithm works only with unique synteny blocks, we filter
all repetitive ones before building the breakpoint graph. Therefore, some
target sequences (generally, short and repetitive contigs) will be ignored.
Some portion of them is put back on the refinement step of the algorithm.

To incorporate these repetitive fragments into final scaffolds you can
turn on new experimental algorithm, which tries to separate different instances
of each single repeat based on their "context" (--repeats option). Depending
on the input assembly, you may get a significant increase in the number of
incorporated fragments, therefore decreasing scaffolds gaps. However,
if there are copy number variations between reference and target genomes,
the algorithm could make some false insertions.


Refinement with Assembly Graph
------------------------------

Ragout uses assembly (overlap) graph to incorporate very short / repetitive
contigs into assembly. First, this graph is reconstructed by overlapping
input contigs/scaffolds (see Sequence Data section for requirements).
Then Ragout scaffolds which are already available are being "threaded"
through this graph to find the true "genome path".

This procedure:

* Increases number of conitgs in output scaffolds
* Improves estimates of distances between contigs

However, sometimes the assembly graph is not accurate: some adjacencies
between contigs could be missing and on the other hand, there also might
be some false-positive adjacencies. This may lead to some incorrectly inserted
contigs. 

For good assemblies (with reasonable coverage, read length etc.) the fraction
of such errors should be very small or even zero. And even they exist,
they are "local" and do not violate the genome structure (probably, most
of them even will not be detected as missassembles by quality assesment
tools like Quast).

As this step may take a lot of time for assemblies with big number of contigs,
you may skip it by specifying "--no-refine" option.


Links File
----------

Ragout outputs information about generated adjacencies in "*.links" file.
It is organized as a table for each scaffold and includes values described below:

* __contig_1__ : first contig in adjacency
* __contig_2__ : second contig in adjacency
* __gap__ : estimated gap between contigs
* __ref_support__ : reference genomes that support this adjacency
* __~>__ : indicates that this adjacency was generated during the refinement procedure

Gap < 0 means overlap on a corresponding value. "~>" does not
apply for assemblies without refinement.


Useful Scripts
--------------

Scripts are located in "scripts" directory

**verify-order.py:**

Tests the correctness of the inferred contigs order if a "true" reference
is available. First, contigs should be mapped on that reference using 
*nucmer* software:

    nucmer --maxmatch --coords reference contigs

Then run the script with the obtained "coords" file:

    scripts/verify-order.py nucmer_coords ord_file