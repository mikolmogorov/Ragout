Installation instructions for Ragout
====================================

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
Download and install it into your system (binary packages are available
for all popular platforms).

If you prefer local installation, the main executable file should be in 
your system executable path ("PATH" environment variable in UNIX).
Also, you can edit "SIBELIA_BIN" variable in "ragout.py" file.

When ready, you can test your installation on provided ready-to-use
examples (see "USAGE.md")