Installation instructions for Ragout
====================================


Build requirements
------------------
* Python 2.7
* C++ compiler with C++11 support (gcc 4.7+ or proper version of Clang)
* Cmake (for building Sibelia)
* Some standard POSIX utilities, such as *wget* or *tar*


Runtime depenencies
-------------------

* python 2.7
* biopython [http://biopython.org]
* networkx [http://networkx.github.io]
* Sibelia [https://github.com/bioinf/Sibelia]


Building and installing
-----------------------

Ragout is distributed as a Python package. Currently, only Python 2.7
is supported. To build and install Ragout run:

	python2.7 setup.py install

This also will install all necessary python dependencies, if neded.
Otherwise, you can install them manually using *pip* or your OS-specific
package manager.

You also can specify installation prefix:

	python2.7 setup.py install --prefix=/usr/local

Linux users would probably prefer to install to "~/.local". 
In this case you may need to add "~/.local/bin" to your
executable path (to run Ragout from any working directory):

	export PATH=$PATH:~/.local/bin

Ragout requires Sibelia for synteny block decomposition.
To instal it, run:

	scripts/install-deps.py --prefix=your_prefix

Do not forget that "your_prefix/bin" also should be in your PATH.
After installation process you can test your installation by running:

	ragout --help

If it works, you can try Ragout on the provided examples (refer to USAGE.md for this)