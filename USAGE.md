Usage
=====

    python ragout.py [-h] -c config -o output_dir [-s] [-g]
    
Supported arguments:

    -h, --help      show the help message and exit
    -c config       Configuration file
    -o output_dir   Output directory
    -s              Skip Sibelia running step
    -g              Refine with assembly graph


Examples:
=========

You can try Ragout on the provided ready-to-use examples:

    python ragout.py -c examples/E.Coli/ecoli.cfg -o examples/E.Coli/out/
    python ragout.py -c examples/H.Pylori/helicobacter.cfg -o examples/H.Pylori/out/
    python ragout.py -c examples/S.Aureus/aureus.cfg -o examples/S.Aureus/out/
    python ragout.py -c examples/V.Cholerea/cholerea.cfg -o examples/V.Cholerea/out/


Input:
======

Ragout takes as input:

- Reference sequences in "fasta" format
- Target assembly in "fasta" format
- Phylogenetic tree for both references and target assembly in "newick" format
- Set of minimum synteny block sizes (one for each iteration)

All these settings should be described in a single config file.
See the example of such file below.


Configuration file
==================

Here is an example of Ragout configuration file:

    REF col=references/COL.fasta
    REF jkd=references/JKD6008.fasta
    REF rf122=references/RF122.fasta
    REF n315=references/N315.fasta

    TARGET usa=usa300_contigs.fasta

    TREE=(rf122:0.0280919,(((usa:0.0151257,col:0.0127906):0.0132464,jkd:0.0439743):0.00532819,n315:0.0150894):0.0150894);

    BLOCK=5000,500,100

Keywords explanation:

- REF: label of reference sequence and it`s relative path
- TARGET: label of target assembly and it`s relative path
- TREE: phylogenetic tree for both references and target assembly
- BLOCK: set of minimum synteny block sizes (one for each iteration)


Output files
============

After running Ragout, an output directory will contain:

* "scaffolds.ord" with the resulting contigs order
* "scaffolds.fasta" with the scaffold sequences (contigs are separated by 11 Ns)
* "scaffolds_refined.ord" with the refined contigs order (if you used -g option)
* "scaffolds_refined.fasta" with the refined scaffold sequences (if you used -g option)


Useful scripts
==============

Scripts are located in "scripts" directory

test.py:

    python test.py nucmer_coords ord_file

Tests the correctness of the infered contigs order, if a closely related reference
is available. Script takes "nucmer" "coords" file as the first argument,
and "ord" file as second (see above).