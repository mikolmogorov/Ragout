Ragout
======

Version: 2.0

Website: http://fenderglass.github.io/Ragout/

[![Build Status](https://travis-ci.org/fenderglass/Ragout.svg?branch=master
)](https://travis-ci.org/fenderglass/Ragout)

<p align="center">
  <img src="http://fenderglass.github.io/Ragout/images/ragout.png" alt="Ragout logo"/>
</p>

Description
-----------
Ragout (Reference-Assisted Genome Ordering UTility)
is a tool for chromosome assembly using multiple references. Given a set of assembly fragments
(contigs/scaffolds) and one or multiple related references (complete or draft),
it produces a chromosome-scale assembly (as a set of scaffolds).

The approach is based on the analysis of genome rearrangements
(like inversions or chromosomal translocations) between the input genomes
and reconstructing the most parsimonious structure of the target genome.

Ragout now supports both small and large genomes (of mammalian scale and complexity).
The assembly of highly polymorphic genomes is currently limited.

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


Publications
------------
- Kolmogorov et al., "Chromosome assembly of large and complex genomes using multiple references",
bioRxiv preprint, 2016

- Kolmogorov et al., "Ragout: A reference-assisted assembly tool for bacterial genomes",
Bioinformatics, 2014


Contacts
--------
Please report any issues directly to the github issue tracker.
Also, you can send your feedback to fenderglass@gmail.com


Acknowledgments
----------------
The work was partially supported by VP Foundation.

We would like to thank:
- Anna Liosnova (benchmarks and useful suggestions)
- Nikolay Vyahhi (testing and useful suggestions)
- Aleksey Gurevich (testing)


Third-party
-----------
Ragout package includes some third-patry software (see INSTALL.md for details)

* Networkx 1.8 Python library [http://networkx.github.io/]
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
