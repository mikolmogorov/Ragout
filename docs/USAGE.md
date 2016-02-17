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

      --solid-scaffolds     do not break input sequences - disables chimera
                            detection module (default: False)
    
      --overwrite           overwrite results from the previous run (default: False)
    
      --repeats             resolve repetitive input sequences (default: False)
    
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


Initial Assembly Requirements
-----------------------------

Ragout takes initial assembly fragments (contigs/scaffolds) as input. 
We performed our tests with SPAdes, ABySS and SOAPdenovo assemblers, others should
work fine too, if their output satisfy the following conditions:

* Assembly coverage should be sufficient (80-90%+)
* *All* contigs/scaffolds output by assembler should be used (including
  short ones)
* For the better performance of the refinment module, contigs/scaffolds
  that were connected in a graph used by assembler should share
  the same sequence on their ends (usually k-mer or (k-1)-mer). This works
  for the most of assemblers which utilize de Bruijn graphs. Currently,
  for other types of assemblers (such as SGA) the performance of
  the refinement module is limited.



Algorithm Overview
-------------------

Ragout works with genomes represented as sequences of synteny blocks
and firstly uses *Sibelia* or *HAL alignment* for this decomposition. 
This step is usually the most time-consuming.

Next, Ragout constructs a breakpoint graph (which reflects adjacencies 
between synteny blocks in the input genomes) and predicts missing adjacencies 
in the target genome (as it is fragmented into contigs/scaffolds, 
some adjacencies are missing). Then, assembly fragments are being joined 
into scaffolds with respect to the inferred adjacencies. This procedure is 
repeated multiple times with the different synteny blocks scale.

Afterwards, a refinement step is performed. Ragout reconstructs
assembly (overlap) graph from the assembly fragments and uses
this graph to insert very short/repetitive fragments into the 
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

After running Ragout, output directory will contain 
(substitute target with the name of your target genome):

* __target_scaffolds.fasta__: scaffolds sequences
* __target_unplaced.fasta__: unplaced input fragments
* __target_scaffolds.links__: description of the input fragments order (see below)
* __target_scaffolds.agp__: a similar description in NCBI AGP format

There will also be some intermediate files related to the run.


Recipe File
-----------

Recipe file describes Ragout run configuration.
Here is an example of such file (a full version, some parameters could be ommited):

    .references = rf122,col,jkd,n315
    .target = usa

    col.fasta = references/COL.fasta
    jkd.fasta = references/JKD6008.fasta
    rf122.fasta = references/RF122.fasta
    n315.fasta = references/N315.fasta
    usa.fasta = usa300_contigs.fasta

    .tree = (rf122:0.02,(((usa:0.01,col:0.01):0.01,jkd:0.04):0.005,n315:0.01):0.01);
    .blocks = small
    .naming_ref = rf122
    

or, if using *HAL* as input, tree, blocks scale and naming reference are inferred automatically

    .references = miranda,simulans,melanogaster
    .target = yakuba
    .hal = genomes/alignment.hal
   

###Parameters description:

Each parameter could be "global" (related to the run) or "local" (for a particular genome).
Global parameters start from dot:

    .global_param_name = value

To set local parameter, use:

    genome_name.param_name = value

###Global parameters

* __references__: comma-separated list of reference names [*required*]
* __target__: target genome name [*required*]
* __tree__: phylogenetic tree in NEWICK format
* __blocks__: synteny blocks scale
* __hal__: path to the alignment in *HAL* format
* __naming_ref__: referece to use for output scaffolds naming

###Local parameters

* __fasta__: path to *FASTA* [default = not set]
* __draft__: indicates that reference is in a draft form (not chromosomes) [default = false]

###Default values

You can change default values for local parameters by assigning the 
parameter value to the special "star" object.
For instance, if all input references except one are in draft form, you can write:

    *.draft = true
    complete_ref.draft = false

###Quick comments

Paths to *FASTA*/*HAL* can be both relative and absolute. 

If you use *Sibelia* for synteny blocks decomposition you must specify 
FASTA for each input genome. If you use *HAL*, everything is taken from it.

Running with Sibelia requires all sequence headers (">gi...") 
among ALL *FASTA* files to be unique.

If you do not specify phylogenetic tree or synteny block scale, 
they will be inferred automatically.


Parameters Description
----------------------

### Phylogenetic tree

Ragout algorithm requires a phylogenetic tree as input. This tree
could be inferred automatically from the breakpoint configuration
of the input genomes. The automatic inference
generally produces a good approximation of a real phylogeny
and is therefore recommended for most of the purposes.
However, if you already have the tree structure from a different source,
you may guide the algorithm with it by setting the corresponding parameter.


### Synteny block scale

Because the decomposition procedure is parameter-dependent, the assembly
is performed in multiple iterations with different synteny block
scale. Intuitively, the algorithm firstly considers only fragments
that are long enough and then insert shorter ones into final scaffolds.

There are two pre-defined scales: "small" and "large". We recommend
"small" for relatively small genomes (bacterial) and large otherwise 
(mammalian). If the parameter is not set, it is automatically inferred
based on input genomes size (recommended).


### Reference genome in draft form

Ragout can use an incomplete assembly (contigs/scaffolds) as a reference.
In such a case you should specify that the reference is in draft from by
setting the corresponding parameter in the recipe file.


### Naming reference

Output scaffolds will be named according to the homology to one single 
reference (naming reference). This reference could be set with a corresponding
recipe parameter, otherwise it would be chosen as the closest reference in the
phylogenetic tree. The naming pattern is as follows. If a scaffold is homologous
to a single reference chromosome "A", it would be named as "chr_A". If there
are multuple homologous chromosomes, for example "A" and "B" (in case of 
chromosomal fusion), it will be named "chr_A_B". If there are multiple
scaffodlds with a same name, the longest one would be chosen as primary,
others will get an extra "_unlocalized" suffix.


Synteny backends
----------------

Ragout has two different options for synteny block decomposition:

* Decomposition with *Sibelia*
* Use of whole genome alignment (in *HAL* format)

You can choose between backends by specifying --synteny (-s) option.
If you use *Sibelia*, you should specify separate *FASTA* file for 
each input genome, while if you work with *HAL*, it is not necessary.

### Sibelia

"Sibelia" is the default option and recommended for bacterial genomes.

### Whole genome alignment in *HAL* format

Alternatively, Ragout can use *HAL* whole genome alignment for synteny blocks
decomposition. This option is recommended for large (over 100MB) genomes, which
*Sibelia* can not process. This alignment could be done by Progressive Cactus 
aligner [https://github.com/glennhickey/progressiveCactus].

The pipeline is as follows: first, align references and target genomes using
Progressive Cactus (you will need phylogenetic tree) and then run Ragout
with the resulted *HAL* file. This would require *HAL Tools* package to be 
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
target sequences (generally, short and repetitive contigs) will be ignored
(some portion of them will be put back during the refinement step of the algorithm).

To incorporate these repetitive fragments into the assembly you can
use the experimental algorithm, which tries to resolve the repetitive contigs
and find their positions in the assmebly (--repeats option). Depending
on nature of the input, you may get a significant increase in the assembly
coverage, therefore decreasing scaffolds gaps. However, if there are copy number 
variations between reference and target genomes, the algorithm could make 
some false insertions.


Chimera detection
-----------------

Ragout can detect chimeric adjacencies inside the input sequences and fix
them by breaking the sequences into parts. The chimera detection algorithm
tries to distinguish such erroneous joins from target-specific adjacencies,
that are not observed in the references. By default, the adjacency which is
not supported by references is considered chimeric, unless there is an
evidence of a rearrangement in the target genome. Sometimes, due to the 
fragmentation of the target genome, such evidence support is missing
for true rearrangements. If you have high quality contigs/scaffolds,
you may choose not to break them at all by specifying --solid-scaffolds option.


Refinement with Assembly Graph
------------------------------

Ragout uses assembly (overlap) graph to incorporate very short / repetitive
contigs into the assembly. First, this graph is reconstructed by overlapping
input contigs/scaffolds (see Sequence Data section for requirements).
Then Ragout scaffolds which are already available are being "threaded"
through this graph to find the true "genome path".

This procedure:

* Increases number of conitgs in output scaffolds
* Improves estimates of distances between contigs

Sometimes assembly graphs are not accurate: some adjacencies
between contigs could be missing and on the other hand, there also might
be some false-positive adjacencies. This may lead to some incorrectly inserted
contigs. However, for good assemblies (with reasonable coverage, read length etc.) 
the fraction of such errors should be minor.

As this step maight be computationaly expensive for big assemblies,
you can skip it by specifying "--no-refine" option.


Links File
----------

Ragout outputs information about generated adjacencies in "*.links" file.
It is organized as a table for each scaffold and includes values described below:

* __sequence__ : input fragment's name and strand (possibly with coordinates in form [start:end])
* __start__ : fragment's position in the scaffold
* __length__ : fragment's length
* __gap__ : gap size between the current and the next fragment
* __support__ : names of the references that support corresponding adjacency

Input fragments are described in a form:

    [+/-]seq_name[start:end]

Sign corresponds to a fragment's strand. If the [start:end] structure is ommited, 
the full fragment is used. A symbol "~>" in support field means the support 
of the assembly graph.


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
