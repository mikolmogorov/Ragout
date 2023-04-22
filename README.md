Ragout
======

Version: 2.3


[![Build Status](https://travis-ci.org/fenderglass/Ragout.svg?branch=master)](https://travis-ci.org/fenderglass/Ragout)

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/ragout.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/ragout)

<p align="center">
  <img src="http://fenderglass.github.io/Ragout/images/ragout.png" alt="Ragout logo"/>
</p>

Description
-----------
Ragout (Reference-Assisted Genome Ordering UTility)
is a tool for chromosome-level scaffolding using multiple references. 
Given initial assembly fragments (contigs/scaffolds) and one or multiple 
related references (complete or draft), it produces a chromosome-scale assembly
(as a set of scaffolds).

The approach is based on the analysis of genome rearrangements
(like inversions or chromosomal translocations) between the input genomes
and reconstructing the most parsimonious structure of the target genome.

Ragout now supports both small and large genomes (of mammalian scale and complexity).
The assembly of highly polymorphic genomes is currently limited.


Manuals
-------

- [Installation instructions](docs/INSTALL.md)
- [Usage](docs/USAGE.md)


Code contributions
------------------
* Mikhail Kolmogorov (St. Petersburg University of the Russian Academy of Sciences, UCSD)
* Pavel Avdeev (St. Petersburg University of the Russian Academy of Sciences)
* Dmitriy Meleshko (St. Petersburg University of the Russian Academy of Sciences)
* Son Pham (UCSD)
* Tatiana Malygina


Publications
------------
* Kolmogorov M, Armstrong J, Raney BJ, Streeter I, Dunn M, Yang F, Odom D, Flicek P, Keane TM, 
Thybert D, Paten B. and Pham S. "Chromosome assembly of large and complex genomes using multiple references" 
Genome research. 2018 [doi:10.1101/gr.236273.118](https://genome.cshlp.org/content/28/11/1720.short)

* Kolmogorov, M., Raney, B., Paten, B. and Pham, S. 
"Ragout—a reference-assisted assembly tool for bacterial genomes" 
Bioinformatics, 2014 [doi:10.1093/bioinformatics/btu280](https://doi.org/10.1093/bioinformatics/btu280)

Contacts
--------
Please report any issues directly to the github issue tracker.
Also, you can send your feedback to mikolmogorov@gmail.com


Acknowledgments
---------------
The work was partially supported by VP Foundation.

We also would like to thank:
* Anna Liosnova (benchmarks and useful suggestions)
* Nikolay Vyahhi (testing and useful suggestions)
* Aleksey Gurevich (testing)


Third-party
-----------
Ragout is using some third-patry software (see INSTALL.md for details):

* Networkx Python library [http://networkx.github.io/]
* Newick parser by Thomas Mailund [https://www.mailund.dk/]
* Sibelia [http://github.com/bioinf/Sibelia]
* HAL Tools [https://github.com/ComparativeGenomicsToolkit/hal]


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
https://zenodo.org/record/2633314/files/ismb_171_supplementary.zip
