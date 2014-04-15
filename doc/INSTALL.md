Installation instructions for Ragout
====================================


Build requirements
------------------
* Python 2.7
* C++ compiler with C++11 support (GCC 4.7+ or proper version of Clang)
* Cmake (for building Sibelia)
* Some standard POSIX utilities, such as *wget* or *tar*


Runtime depenencies
-------------------

* Python 2.7
* biopython [http://biopython.org]
* networkx [http://networkx.github.io]
* Sibelia [https://github.com/bioinf/Sibelia]


You can install Ragout as a Python package, which is recommended.
Alternatively, you can build and run it from a source directory.


Package installation (recommended)
----------------------------------

To install Ragout as a Python package

	python2.7 setup.py install

If you don't have permission to install software on your system, you can 
install into another directory using the --user, --prefix, or --home flags to setup.py.

	python2.7 setup.py install --user
	or
	python2.7 setup.py install --prefix=~/.local
	or
	python2.7 setup.py install --home=~

If you didn't install in the standard Python site-packages directory you will 
need to set your PYTHONPATH variable to the alternate location. 
See http://docs.python.org/2/install/index.html#search-path for further details.

After installation with custom prefix you may need to add the corresponding 
"bin" directory to your executable path (to run Ragout from any working directory). 
For example, if your prefix was "~/.local", run:

	export PATH=$PATH:~/.local/bin

setup.py script also will install all necessary Python dependencies, if neded.
Alternatively, you can install them manually using *pip* or your OS-specific
package manager:

	pip install biopython networkx
	or
	sudo apt-get install biopython networkx

After installation process you can test your installation by running:

	ragout --help

If it works, you can try Ragout on the provided examples (refer to USAGE.md for this)


Using without installation
--------------------------

To build Ragout in source directory, run:

	python2.7 setup.py inplace

This will build all necessary modules and create "ragout_local" executable.
In this case, you should manually install all dependencies using *pip*
or your OS-specific package manager as it is written below.


Sibelia
-------

Ragout requires Sibelia for synteny block decomposition.
To instal it, run:

	scripts/install-deps.py --prefix=your_prefix

Do not forget that "your_prefix/bin" also should be in your PATH.