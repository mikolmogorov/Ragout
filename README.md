Ragout
======

Version: 0.2 beta
Release date: 15 Apr 2014

Website: http://fenderglass.github.io/Ragout/



       	         (
		     )    )
		  _.(--"("""--.._
		 /, _..-----).._,\
		|  `'''-----'''`  |
		 \               /
		  '.           .'
		    '--.....--'

Description
-----------
Ragout (Reference-Assisted Genome Ordering UTility) is a tool for
assisted assembly using multiple references. It takes a short read
assembly (a set of contigs), a set of related references
and a corresponding phylogenetic tree and then assembles the contigs into
scaffolds.

The benefits of assembly with multiple references become significant,
when those references have structural variations compared to the target
genome. Even if each reference is structurally divergent, it is possible
to assemble the target into the correct set of scaffolds. Enlarge your
contigs with Ragout!

The current version of Ragout is limited to genomes of bacterial size,
but we are working on expanding it to mammalian-scaled ones.


Install
-------
See *docs/INSTALL.md* file.

Usage
-----
See *docs/USAGE.md* file.


Authors
-------
- Mikhail Kolmogorov (St. Petersburg University of the Russian Academy of Sciences)
- Pavel Avdeev (St. Petersburg University of the Russian Academy of Sciences)
- Dmitriy Meleshko (St. Petersburg University of the Russian Academy of Sciences)
- Son Pham (University of California, San Diego)


Citation
--------
- Mikhail Kolmogorov, Brian Raney, Benedict Paten, and Son Pham. 
"Ragout: A reference-assisted assembly tool for bacterial genomes",
Bioinformatics, 2014


ISMB 2014 supplementary
-----------------------

Supplementary materials for ISMB submission could be found at:
https://drive.google.com/file/d/0B1pUguR1yn7TMjNpX09JdFphT3c/edit?usp=sharing


Contacts
--------
Please report any bugs directly to the issue tracker of this project.
Also, you can send your feedback at fenderglass@gmail.com


Acknowledgements
----------------
The work was supported by VP Foundation.

We would like to thank:
- Nikolay Vyahhi (testing and some useful suggestions)
- Aleksey Gurevich (testing)


Third-party
-----------
Ragout uses some third-patry software;

* Networkx Python library [http://networkx.github.io/]
* Newick Python parser [http://www.daimi.au.dk/~mailund/newick.html]
* Sibelia [http://github.com/bioinf/Sibelia]


License
-------
Ragout is a distribution package that contains several
other pieces of free software. Most of this
software is completely free to redistribute, but some such as Sibelia
or Newick parser are released under the GPL. We therefore release
Ragout distribution under the GPL and note that the
licenses of the constituent packages can be viewed in their
subfolders. (see *LICENSE* file)