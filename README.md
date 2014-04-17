Ragout
======

Version: 0.2 beta
Release date: 15 Apr 2014

Github: https://github.com/fenderglass/Ragout
PyPI: https://pypi.python.org/pypi/ragout


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
- Mikhail Kolmogorov (St. Petersburg University of the Russian Academy of
Sciences)
- Son Pham (University of California, San Diego)


Citation
--------
- Mikhail Kolmogorov, Brian Raney, Benedict Paten, and Son Pham. 
"Ragout: A reference-assisted assemble tool for bacterial genomes" 
(accepted to ISMB 2014)


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


Licence
-------
The program is distributed under GNU GPL v2 license.