Installation instructions for Ragout
====================================

Ragout is written in python and requires version 2.7 to run.
Also, there are third-party dependencies:

* biopython [http://biopython.org]
* networkx [http://networkx.github.io]

You can install them with your OS-specific package manager,
e.g. in Ubuntu:

$ sudo apt-get install python-biopython python-networkx

or using pip package manager:

$ pip install biopython networkx

Ragout can use Sibelia for genomes decomposition on synteny blocks.
"Sibelia" executable should be in your "PATH" environment variable,
otherwise you can edit "SIBELIA_BIN" variable in "ragout.py" file. 

If you use Python 2.5 or 2.6, you need to install additional Python
packages:

* argparse [https://pypi.python.org/pypi/argparse]
