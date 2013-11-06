Usage
=====

    python ragout.py [-h] -c config -o output_dir [-s]
    
Supported arguments:

    -h, --help     show the help message and exit
    -c config      Configuration file
    -o output_dir  Output directory
    -s             Skip Sibelia running step

You can try Ragout on provided ready-to-use examples:

    python ragout.py -c examples/E.Coli/ecoli.cfg -o examples/E.Coli/out/
    python ragout.py -c examples/H.Pylori/helicobacter.cfg -o examples/H.Pylori/out/
    python ragout.py -c examples/S.Aureus/aureus.cfg -o examples/S.Aureus/out/

Configuration file
==================

Unless you want to use only our examples, Ragout should be provided with meaningful config file. Here is an example of such file:

    REF col=references/COL.fasta
    REF jkd=references/JKD6008.fasta
    REF rf122=references/RF122.fasta
    REF n315=references/N315.fasta
    
    TARGET usa=usa300_contigs.fasta
    
    TREE=(jkd:0.0359,(((col:0,target:0.00028)0.91:0.00153,n315:0.00323)1:0.03131,rf122:0.0277):0.00135);
    
    BLOCK=5000,500,100

Supported keywords:

- REF indicates reference sequence (in fasta format).
- TARGET refers to file containing contig sequences (in fasta format).
- TREE represents phylogenetic tree of references and target (in Newick format).
- BLOCK corrseponds to set of minimum synteny block sizes for each iteration.

Output files
============

After running Ragout, an output directory will contain:

- scaffolds.fasta with resulting assembly.
- ...

Useful scripts
==============

Scripts are located in "scripts" directory

    test.py nucmer_coords ord_file 

Tests correctness of inferred contigs order, if closely related reference
sequence available. It takes nucmer`s coords file as the first argument,
and ord file, which was produced by Ragout.
