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
    
      --refine              enable refinement with assembly graph (default:
                            False)

      --solid-scaffolds     do not break input sequences - disables chimera
                            detection module (default: False)
    
      --overwrite           overwrite results from the previous run (default: False)
    
      --repeats             enable repeat resolution algorithm (default: False)
    
      --debug               enable debug output (default: False)
    
      -t THREADS, --threads THREADS
                            number of threads for synteny backend (default: 1)
    
      --version             show program's version number and exit


Examples
--------

You can try Ragout on the provided ready-to-use examples:

    ./ragout.py examples/E.Coli/ecoli.rcp --outdir examples/E.Coli/out/ --refine
    ./ragout.py examples/H.Pylori/helicobacter.rcp --outdir examples/H.Pylori/out/ --refine
    ./ragout.py examples/S.Aureus/aureus.rcp --outdir examples/S.Aureus/out/ --refine
    ./ragout.py examples/V.Cholerae/cholerae.rcp --outdir examples/V.Cholerae/out/ --refine


Algorithm Overview
-------------------

Ragout first uses *Sibelia* or *HAL alignment* for to decompose the input
genomes into the sequences of synteny blocks -- this step is usually 
the most time-consuming.

Using the synteny information, Ragout infers the phylogenetic tree
of the input genomes (if it is not given as input).
Next, Ragout constructs the breakpoint graph (which reflects adjacencies 
between the synteny blocks in the input genomes) and recovers the missing adjacencies 
in the target genome (as it is fragmented, some adjacencies are missing). 
Then, assembly fragments are joined into scaffolds. The final
chromosomes are constructed as a consensus of scaffolds built
with different synteny block scales.

Finally, an optional refinement step is performed. Ragout reconstructs
assembly (overlap) graph from the assembly fragments and uses
this graph to insert very short/repetitive fragments into the 
assembly.


Input
-----

Ragout takes as input:

* Reference genomes [in *FASTA* format or packed into *HAL*]
* Target assembly in [in *FASTA* format or packed into *HAL*]

Optionally, you can add:

* Phylogenetic tree with the reference and target genomes [in *NEWICK* format]
* Synteny block scale

All these parameters should be described in a single recipe file (see below)


Output
------

After running Ragout, output directory will contain:

* __target_scaffolds.fasta__: scaffolds
* __target_unplaced.fasta__: unplaced input sequences
* __target_scaffolds.links__: the order and orientation of the input sequences in scaffolds (see below)
* __target_scaffolds.agp__: same as below, but in NCBI AGP format


Recipe File
-----------

A recipe file describes the Ragout run configuration.
Here is an example of such file (full version, some parameters could be ommited):

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
* __naming_ref__: reference to use for output scaffolds naming

###Local parameters

* __fasta__: path to *FASTA* [default = not set]
* __draft__: indicates that reference is in a draft form (not chromosomes) [default = false]

###Default values

You can change default values of the local parameters by assigning the 
parameter value to the special "star" object:
for instance, if all input references except one are in a draft form, you can write:

    *.draft = true
    complete_ref.draft = false

###Quick comments

Paths to *FASTA*/*HAL* can be both relative and absolute. 

If you use *Sibelia* for synteny blocks decomposition you must specify 
FASTA for each input genome. If you use *HAL*, all sequences will be taken from it.

*Sibelia* requires all sequence headers (">gi...") 
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

Output scaffolds will be named according to a homology to one of the input 
references (naming reference). This reference can be set with the corresponding
recipe parameter, otherwise it will be chosen as the closest reference in the
phylogenetic tree. The naming pattern is as follows. If a scaffold is homologous
to a single reference chromosome "A", it will be named as "chr_A". If there
are multiple homologous chromosomes, for example "A" and "B" (in case of 
chromosomal fusion), it will be named "chr_A_B". If there are multiple
scaffolds with a same name, the longest one would be chosen as primary,
others will get an extra "_unlocalized" suffix.


Synteny Backends
----------------

Ragout has two different options for synteny block decomposition:

* Decomposition with *Sibelia*
* HAL alignment produced by Progressive Cactus

You can choose between backends by specifying --synteny (-s) option.


### Sibelia

"Sibelia" is the default option and recommended for bacterial genomes.

### Whole genome alignment in *HAL* format

Alternatively, Ragout can use *HAL* whole genome alignment for synteny blocks
decomposition. This option is recommended for large (over 100MB) genomes, which
*Sibelia* can not process. This alignment is done by Progressive Cactus 
aligner [https://github.com/glennhickey/progressiveCactus]. Currently, we do not
provide bindings for running Progressive Cactus from Ragout, as the procedure
might vary for different setups. Ragout starts with the alignment
result in *HAL* format, *HAL tools* should be installed in your system.


### MAF backend is deprecated

The support of MAF synteny backend is deprecated, because it is more 
convenient to work directly with *HAL*, which is a default output of
Progressive Cactus.


Repeat Resolution
-----------------

As the main Ragout algorithm works only with unique synteny blocks, we filter
all repetitive blocks before building the breakpoint graph. Therefore, some
target sequences (generally, short and repetitive contigs) will be ignored
(some of them could be put back during the refinement step below).

To incorporate these repetitive fragments into the assembly, you can
use the optional algorithm, which tries to resolve the repetitive contigs
and find their positions in the assembly ('--repeats' option). Depending
on the dataset, you may get a significant increase in the assembly
coverage, therefore decreasing scaffolds gaps. However, if there are copy number 
variations between the reference and target genomes, the algorithm could make 
some false insertions.


Chimera Detection
-----------------

Ragout detects chimeric adjacencies inside the input sequences and fixes
them by breaking the sequences into parts. The chimera detection algorithm
tries to distinguish such erroneous joins from target-specific adjacencies,
that are not observed in the references. By default, the adjacency which is
not supported by references is considered chimeric, unless there is an
evidence of a rearrangement in the target genome. Sometimes, due to the 
fragmentation of the target genome, the evidence support is missing.
If you have high quality contigs/scaffolds,
you may choose to turn chimera detection off by specifying '--solid-scaffolds' option.


Refinement with the Assembly Graph
----------------------------------

Ragout optionally uses assembly (overlap) graph to incorporate very short / repetitive
contigs into the assembly ('--refine' option). First, this graph is reconstructed by overlapping
input contigs/scaffolds. Then current Ragout scaffolds are "threaded"
through this graph to find the true "genome path".
This procedure increases number of contigs in output scaffolds and also
improves the scaffold gaps estimates. Sometimes assembly graphs are not very accurate, 
which may lead to incorrectly inserted contigs. However, for the most bacterial 
assemblies the fraction of errors should be minor. This procedure is generally recommended
for bacterial assemblies, however, the effect is usually minor for large genomes
because of complications with assembly graph reconstruction.


Links File
----------

Ragout outputs information about generated adjacencies in "*.links" file.
It is organized as a table for each scaffold with the values below:

* __sequence__ : input fragment's name and strand (possibly with coordinates in form [start:end])
* __start__ : fragment's position in the scaffold
* __length__ : fragment's length
* __gap__ : gap size between the current and the next fragment
* __support__ : reference support of the corresponding adjacency

Input fragments are described in a form:

    [+/-]seq_name[start:end]

The sign corresponds to the fragment's strand. The [start:end] structure is omitted
if the full fragment is used. A symbol "~>" in support field corresponds to the 
assembly graph support 


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
