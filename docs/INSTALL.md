Installation instructions for Ragout
====================================


Availability
------------
Ragout is tested under Mac OS and Linux. While it *should* work
under Windows, we currently do not provide an official support.


Build requirements
------------------
* Python 2.7 (with developer headers)
* C++ compiler with C++0x support (GCC 4.6+ / Clang 3.2+ / Apple Clang 4.2+)
* GNU make / Cmake (for building Sibelia)


Runtime depenencies
-------------------

* Python 2.7
* biopython [http://biopython.org]
* networkx 1.8+ [http://networkx.github.io]
* Sibelia [https://github.com/bioinf/Sibelia]


You can install Ragout as a Python package, which is recommended.
Alternatively, you can build and run it from a source directory.

Buildng Ragout requires Python headers, which should be
explicitly installed on some OS. For instance, Ubuntu
users should check that package "python-dev" is installed.


Installation from PyPI (recommended)
------------------------------------

The easies way to install Ragout is to use Python package index database:

	pip install ragout

Note, that this may require superuser privileges:

	sudo pip install ragout

If you do not have *pip* installed, you can get it from here:
http://www.pip-installer.org/


Binary packages
---------------

Pre-compiled binary packages for Linux and Mac OS are available at:
https://pypi.python.org/pypi/ragout as binary *eggs*.
You can install them via *easy_install*:
http://pythonhosted.org/setuptools/easy_install.html


Installing from source
----------------------

To install Ragout as a Python package, run:

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
After installation process you can test your installation by running:

	ragout --help

If it works, you can try Ragout on the provided examples (refer to USAGE.md for this)


Using without installation
--------------------------

To build Ragout inside the source directory, run:

	python2.7 setup.py inplace

In this case, you should manually install all dependencies using *pip*
or your OS-specific package manager:

	pip install biopython networkx
	or
	sudo apt-get install biopython networkx


Sibelia
-------

Ragout requires Sibelia for synteny block decomposition.
You can download and install it from the website: https://github.com/bioinf/Sibelia

Otherwise, you can use our script for a quick installation:

	sudo scripts/install-sibelia.py
	or
	scripts/install-sibelia.py --prefix=your_prefix

Alternatively, if you have installed Ragout with *pip* and do not have
"scripts" directory:

	curl https://raw.github.com/fenderglass/Ragout/master/scripts/install-sibelia.py -o - | python

	or, if you want more control:

	curl https://raw.github.com/fenderglass/Ragout/master/scripts/install-sibelia.py -o - | [sudo] python [- --prefix=your_prefix]

Do not forget that "your_prefix/bin" folder also should be in your PATH.
Alternatively, you can set SIBELIA_INSTALL variable to directory
containing *Sibelia* excecutable.


Troubleshooting
---------------

Q: I get compilation error "Python.h: No such file or directory":

A: You do not have Python developer headers installed. On some
systems/distributions you have to install them explicitly, i.e. in Ubuntu:
	
	sudo apt-get install python-dev


Q: Multiple errors during compilation, possibly with 
"unrecognized command line option '-std=c++0x'" message:

A: Probably your compiler is too old and does not support C++0x. Minimum
versions of GCC and Clang are mentioned above.

Q: clang error: unknown argument: '-mno-fused-madd'

A: This is a python/clang issue on Mac OS. To resolve it, update
your Python to version 2.7.6+ 