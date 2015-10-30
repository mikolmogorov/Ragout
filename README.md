Ragout
======

Version: 1.2

Release date: 28 Sep 2015

Website: http://fenderglass.github.io/Ragout/

[![Build Status](https://travis-ci.org/fenderglass/Ragout.svg?branch=master
)](https://travis-ci.org/fenderglass/Ragout)



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
Ragout (Reference-Assisted Genome Ordering UTility)
is a tool for reference-assisted assembly. Given a set of initial sequences 
(contigs/scaffolds) and one or multiple references (complete or draft) as input
it produces a chromosome-scaled assembly (as a set of gapped scaffolds).

The approach is based on the analysis of medium- and large-scale rearrangements
(like inversions or chromosomal translocations) between the input genomes
using breakpoint graph.

The first version of Ragout was limited to bacterial genomes only,
but currently there is an experimental support of mammalian-scaled genomes as well.


Install
-------
See *docs/INSTALL.md* file.

Usage
-----
See *docs/USAGE.md* file.


Authors
-------
- Mikhail Kolmogorov (St. Petersburg University of the Russian Academy of Sciences, UCSD)
- Pavel Avdeev (St. Petersburg University of the Russian Academy of Sciences)
- Dmitriy Meleshko (St. Petersburg University of the Russian Academy of Sciences)
- Son Pham (UCSD)


Citation
--------
- Mikhail Kolmogorov, Brian Raney, Benedict Paten, and Son Pham. 
"Ragout: A reference-assisted assembly tool for bacterial genomes",
Bioinformatics, 2014


Contacts
--------
Please report any problems directly to the github issue tracker.
Also, you can send your feedback to fenderglass@gmail.com


Acknowledgements
----------------
The work was partially supported by VP Foundation.

We would like to thank:
- Anna Liosnova (benchmarks and useful suggestions)
- Nikolay Vyahhi (testing and useful suggestions)
- Aleksey Gurevich (testing)


Third-party
-----------
Ragout package includes some third-patry software (see INSTALL.md for details)

* Networkx Python library [http://networkx.github.io/]
* Newick Python parser [http://www.daimi.au.dk/~mailund/newick.html]
* Sibelia [http://github.com/bioinf/Sibelia]


License
-------
Ragout itself is distributed under BSD license, but the package also contains
some third-party software. Most of this software is completely free to redistribute,
but some such as Sibelia or Newick parser are released under the GPL. We therefore release
Ragout distribution under the GPL and note that the licenses of the constituent
packages can be viewed in their subfolders. (see *LICENSE* file)


ISMB 2014 supplementary
-----------------------

Supplementary materials for ISMB submission could be found at:
https://drive.google.com/file/d/0B1pUguR1yn7TMjNpX09JdFphT3c/edit?usp=sharing
