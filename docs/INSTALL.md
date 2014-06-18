Installation instructions for Ragout
====================================


Availability
------------
Ragout is tested under Mac OS and Linux. While it *should* work
under Windows, we currently do not provide an official support.


Build requirements
------------------
* Python 2.7 (with developer headers installed)
* C++ compiler with C++0x support (GCC 4.6+ / Clang 3.2+ / Apple Clang 4.2+)
* GNU make / Cmake (for building Sibelia)


Runtime depenencies
-------------------

* Python 2.7
* biopython [http://biopython.org]
* networkx 1.8+ [http://networkx.github.io]
* Sibelia [https://github.com/bioinf/Sibelia] or Progresssive Cactus [https://github.com/glennhickey/progressiveCactus]


Buildng Ragout requires Python headers, which should be
explicitly installed on some OS. For instance, Ubuntu
users should check that package "python-dev" is installed.

You have two options for the installation: First, you can build Ragout
and run it from the distribution directory without installation. 
Alternatively, you can install it as a Python package to integrate into
your system. The second option is recommended, however, it may require
superuser privilegies.


Using from distribution folder
------------------------------

To build Ragout inside the distribution directory, run:

	python2.7 setup.py inplace

With this type of setup, you should manually install all dependencies using *pip*
or your OS-specific package manager:

	pip install biopython networkx

or

	sudo apt-get install biopython networkx

After the installation process you can test your installation by running:

	./ragout-local --help

If it works, you can try Ragout on the provided examples (refer to USAGE.md for this)


Installing as a package
-----------------------

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

After the installation process you can test your installation by running:

	ragout --help

If it works, you can try Ragout on the provided examples (refer to USAGE.md for this)


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


progressiveCactus
-----------------

First, download and build Progressive Cactus: https://github.com/glennhickey/progressiveCactus
Then set "CACTUS_INSTALL" variable pointing to cactus installation directory:


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