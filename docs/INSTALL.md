Installation Instructions for Ragout
====================================

Availability
------------
Ragout is available for Mac OS and Linux.


Bioconda distribitions
----------------------

Ragout binary releases are available through Bioconda:

    conda install -c bioconda ragout


Build Requirements
------------------
* C++ compiler with C++0x support (GCC 4.6+ / Clang 3.2+ / Apple Clang 4.2+)
* GNU make 
* Cmake


Runtime Depenencies
-------------------

* Python (2.7 or 3.5+)
* Sibelia [http://github.com/bioinf/Sibelia]
* python-networkx == 2.2
* HAL Tools (optionally) [https://github.com/ComparativeGenomicsToolkit/hal]


Local installation
------------------

If you don't want to use bioconda release, you can build
Ragout repository clone and run it locally without installing
into system. To do this, perform:

    git clone https://github.com/fenderglass/Ragout.git
	cd Ragout
	python setup.py build
	pip install -r requirements.txt --user
    python scripts/install-sibelia.py

This will also build and install Sibelia and all Python dependencies.
See below for HAL installation instructions.

Once installed, you can invoke Ragout from the cloned directory by using:

    bin/ragout

System installation
-------------------

To integrate Ragout into your system, run:

    git clone https://github.com/fenderglass/Ragout.git
	cd Ragout
	python setup.py build
    python setup.py install

This assumes that you already have python-networkx package
installed into your system (using the respective package manager).
Sibelia / HAL tools should also be installed / integrated separately.


HAL Tools
---------

HAL alignment produced by Cactus could be used for synteny 
blocks decomposition instead of Sibelia (recommended for large genomes). 

If you want to use HAL alignment as input,
you need to install HAL Tools package [https://github.com/ComparativeGenomicsToolkit/hal]
as it is described in the manual. Do not forget to properly set PATH and PYTHONPATH
environment variables.


Troubleshooting
---------------

Q: Many compilation errors, possibly with 
"unrecognized command line option '-std=c++0x'" message:

A: Probably your compiler is too old and does not support C++0x. Make
sure you have at least GCC 4.6+ / Clang 3.2+


Q: "libstdc++.so.6: version `CXXABI_1.3.5' not found" or similar error when running

A: Ensure that the version of libc++ that was used to compile Ragout is similar
to one the you currently using. You can specify an extra search path
to a specific library by setting "LD_LIBRARY_PATH" variable.
