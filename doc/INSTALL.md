Installation instructions for Ragout
====================================

Ragout is written in Python and does not require any preparations.
However, there are some third-party dependencies described below
which need to be installed.

Quick dependencies
------------------

* python 2.7
* biopython [http://biopython.org]
* networkx [http://networkx.github.io]
* Sibelia [https://github.com/bioinf/Sibelia]

Python
------

Ragout is written in python and requires version 2.7 to run.
Also, there are third-party dependencies:

* biopython [http://biopython.org]
* networkx [http://networkx.github.io]

You can install them with your OS-specific package manager,
e.g. in Ubuntu:

	$ sudo apt-get install python-biopython python-networkx

or using pip package manager:

	$ pip install biopython networkx

If you use Python 2.5 or 2.6, you need to install additional Python
packages:

* argparse [https://pypi.python.org/pypi/argparse]

Sibelia
-------

Ragout uses *Sibelia* for decomposing genomes into synteny blocks.
You may install it by running:

	$ python bin/install-deps.py

Building process requires *Cmake* as well as some standard UNIX
executables like *wget* or *tar*.
