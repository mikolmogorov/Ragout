Usage Instructions for Ragout
=============================

Quick Usage
-----------

    usage: ragout [-h] [-o output_dir] [-s {sibelia,maf,hal}] [--refine]
                  [--solid-scaffolds] [--overwrite] [--repeats] [--debug]
                  [-t THREADS] [--version]
                  recipe_file
    
    Chromosome assembly with multiple references
    
    positional arguments:
      recipe_file           path to recipe file
    
    optional arguments:
      -h, --help            show this help message and exit
      -o output_dir, --outdir output_dir
                            output directory (default: ragout-out)
      -s {sibelia,maf,hal}, --synteny {sibelia,maf,hal}
                            backend for synteny block decomposition (default:
                            sibelia)
      --refine              enable refinement with assembly graph (default: False)
      --solid-scaffolds     do not break input sequences - disables chimera
                            detection module (default: False)
      --overwrite           overwrite results from the previous run (default:
                            False)
      --repeats             enable repeat resolution algorithm (default: False)
      --debug               enable debug output (default: False)
      -t THREADS, --threads THREADS
                            number of threads for synteny backend (default: 1)
      --version             show program's version number and exit


Examples
--------

You can try Ragout on the provided ready-to-use examples:

    bin/ragout examples/E.Coli/ecoli.rcp --outdir examples/E.Coli/out/ --refine
    bin/ragout examples/H.Pylori/helicobacter.rcp --outdir examples/H.Pylori/out/ --refine
    bin/ragout examples/S.Aureus/aureus.rcp --outdir examples/S.Aureus/out/ --refine
    bin/ragout examples/V.Cholerae/cholerae.rcp --outdir examples/V.Cholerae/out/ --refine


Algorithm Overview
-------------------

Ragout first uses Sibelia or HAL/MAF alignment for to decompose the input
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

Ragout needs as input:

* Reference genomes [in FASTA format]
* Target assembly in [in FASTA format]
* Optionally, a phylogenetic tree with the reference and target genomes

Alternatively, to process larger genomes (>100Mb) you will need to use 
HAL/MAF alignment produced by Cactus.

All input parameters are be described in a single configuration file (see below)


Output
------

After running Ragout, the output directory will contain (where "target" is the name
of your target genome).

* __target_scaffolds.fasta__: assembled scaffolds
* __target_unplaced.fasta__: unplaced input fragments
* __target_scaffolds.links__: the order and orientation of the input fragments in assembled scaffolds
* __target_scaffolds.agp__: same as above, but in NCBI AGP format


Configuration (Recipe) File
---------------------------

A recipe file describes the Ragout run configuration.
Here is an explicit example (some parameters are optional):

    #reference and target genome names (required)
    .references = rf123,col,jkd,n315
    .target = usa
    
    #phylogenetic tree for all genomes (optional)
    .tree = (rf122:0.02,(((usa:0.01,col:0.01):0.01,jkd:0.04):0.005,n315:0.01):0.01);
    
    #paths to genome fasta files (required for Sibelia)
    col.fasta = references/COL.fasta
    jkd.fasta = references/JKD6008.fasta
    rf122.fasta = references/RF122.fasta
    n315.fasta = references/N315.fasta
    usa.fasta = usa300_contigs.fasta
    
    #synteny blocks scale (optional)
    .blocks = small
    
    #reference to use for scaffold naming (optional)
    .naming_ref = rf122
    
if using HAL as input:

    .references = miranda,simulans,melanogaster
    .target = yakuba
    
    #HAL alignment input. Sequences will be extracted from the alignment
    .hal = genomes/alignment.hal

or, using MAF as input:
   
    .references = miranda,simulans,melanogaster
    .target = yakuba
    
    .maf = alignment.maf
    miranda.fasta = references/miranda.fasta
    simulans.fasta = references/simulans.fasta
    melanogaster.fasta = references/melanogaster.fasta
    yakuba.fasta = yakuba.fasta

Each configuration parameter could be "global" (related to the run) or 
"genomic" (for a particular genome). Global parameters start from dot:

    .global_param_name = value

To set a genomic parameter, use:

    genome_name.param_name = value

### Global parameters

* __references__: comma-separated list of reference names [required]
* __target__: target genome name [required]
* __tree__: phylogenetic tree in NEWICK format
* __blocks__: synteny blocks scale
* __hal__: path to the alignment in HAL format
* __maf__: path to the alignment in MAF format
* __naming_ref__: reference to use for output scaffolds naming

If you do not specify phylogenetic tree or synteny block scale, 
they will be inferred automatically.

### Genomic parameters

* __fasta__: path to FASTA [default = not set]
* __draft__: indicates that reference is in a draft form (not chromosomes) [default = false]

Paths to FASTA/HAL/MAF can be absolute or relative to the recipe file. 

If you use Sibelia for synteny blocks decomposition you must specify 
FASTA for each input genome. If you use HAL, sequnces will be extracted 
from the alignment.

Sibelia requires all FASTA sequence identifiers (">gi...") within ALL files to be unique.

You can use wildcards to set the genomic parameters. 
For instance, if all input references except one are in a draft form:

    *.draft = true
    complete_ref.draft = false


Parameters Description
----------------------

### Phylogenetic tree

Ragout algorithm requires a phylogenetic tree as input. If the tree
if not provided, if will be inferred automatically from the 
breakpoint configuration of the input genomes. The automatic inference
generally produces a good approximation of a real phylogeny
and is therefore recommended for the most runs.


### Synteny block scale

The assembly is performed in multiple iterations with different synteny block
scales. Intuitively, the algorithm initially considers only long and reliable
synteny blocks and then use the shorter ones to fill gaps in final scaffolds.

There are two pre-defined block sets: "small" and "large". We recommend
"small" for relatively short genomes (bacterial) and "large" otherwise 
(>100Mb). If the parameter is not set, it is automatically inferred
based on the sizes of input genomes. You may also use
custom set of block sizes, for example:

    .blocks = 50000,5000


### Reference genome in draft form

Ragout can use an incomplete assembly (contigs/scaffolds) as a reference.
In this case you should set the corresponding parameter in the recipe file
as shown above.


### Naming reference

Output scaffolds will be named according to homology to one of the input 
references ("naming reference"). This reference can be set with the corresponding
recipe parameter, otherwise it will be chosen as the closest reference in the
phylogenetic tree. 

The naming rule is as follows. If a scaffold is homologous
to a single reference chromosome "A", it will be named as "chr_A". If there
are multiple homologous chromosomes, for example "A" and "B" (in case of 
chromosomal fusion), it will be named "chr_A_B". If there are multiple
scaffolds with a same name, the longest one would be chosen as primary,
while all other will get "unlocalized" suffix.


Synteny Block Reconstruction
----------------------------

Ragout has two different options for synteny block decomposition:

* Decomposition with Sibelia
* HAL/MAF alignment produced by Cactus

You can choose between backends by specifying --synteny (-s) option.


### Sibelia

Sibelia is the default option and recommended for bacterial and small eukaryotic genomes.
Ragout automatically runs Sibelia in the beginning.

### Whole genome alignment in HAL/MAF format

Alternatively, Ragout can use HAL whole genome alignment for synteny blocks
decomposition. This option is recommended for larger (over 100Mb) genomes, which
Sibelia can not process. This alignment first should be done using Cactus 
aligner [https://github.com/ComparativeGenomicsToolkit/cactus]. 
Afterwards, run Ragout with the produced HAL file 
("HAL tools" package should be installed in your system).

If you are having troubles with coupling Ragout and HAL tools,
you can manually convert HAL to MAF, and use the MAF alignment
as input for Ragout (using `-s maf` option). In this case,
you should always specify paths to FASTA files for each genome
(as in case of using Sibelia).


Repeat Resolution
-----------------

As the main Ragout algorithm works only with unique synteny blocks, we filter
all repetitive blocks before building the breakpoint graph. Therefore, some
target sequences (generally, short and repetitive contigs) will be ignored
(some of them could be put back during the refinement step below).

To incorporate these repetitive fragments into the final assembly, you can
use the optional repeat resolution algorithm ('--repeats' option). Depending
on the dataset, you may get a significant increase in the assembly
coverage, therefore decreasing scaffolds gaps. However, if there are copy number 
variations between the reference and target genomes, the algorithm could make 
some false insertions.


Chimera Detection
-----------------

Ragout detects chimeric adjacencies inside the input sequences and brakes them. 
By default, an adjacency which is not supported by references is considered chimeric, 
unless there is a clear evidence of a rearrangement in the target genome (from
breakpoint analysis). Sometimes, due to the fragmentation of the target genome, 
the breakpoint support is missing. If you have high quality contigs/scaffolds,
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

Ragout outputs information about generated adjacencies in ".links" file.
It is organized as a table for each scaffold with the values below:

* __sequence__ : input fragment's name, strand and coordinates (see below)
* __start__ : fragment's position in the scaffold
* __length__ : fragment's length
* __gap__ : gap size between the current and the next fragment
* __support__ : reference support for the corresponding adjacency

Input fragments are described in a form:

    [+/-]seq_name[start:end]

The sign corresponds to the fragment's strand. The [start:end] structure is omitted
if the full fragment is used. A symbol "~>" in support field corresponds to
support from the assembly graph.


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
