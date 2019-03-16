Installation Instructions for Ragout
====================================

Availability
------------
Ragout is available for Mac OS and Linux.


Bioconda distribitions
----------------------

Ragout binary releases are available through Bioconda:

    conda install ragout


Build Requirements
------------------
* C++ compiler with C++0x support (GCC 4.6+ / Clang 3.2+ / Apple Clang 4.2+)
* GNU make 
* Cmake


Runtime Depenencies
-------------------

* Python 2.7
* Sibelia [http://github.com/bioinf/Sibelia]
* HAL Tools [https://github.com/glennhickey/hal] (alternatively to Sibelia)


Building
--------

To build Ragout binaries, type:
    
        python setup.py build

You will also need either Sibelia or HAL Tools installed

To build and install Sibelia, use:

        python scripts/install-sibelia.py

If you already have Sibelia installed into your system, it will
be picked up automatically.

Optionally, you may isntall Ragout into your system by typing:

        python setup.py install


HAL Tools
---------

HAL alignment produced by Progressive Cactus could be used for synteny 
blocks decomposition instead of Sibelia (recommended for large genomes). 

If you want to use HAL alignment as input,
you need to install HAL Tools package [https://github.com/glennhickey/hal]
as it is described in the manual. Do not forget to properly set PATH and PYTHONPATH
environment variables.


Troubleshooting
---------------

Q: Many compilation errors, possibly with 
"unrecognized command line option '-std=c++0x'" message:

A: Probably your compiler is too old and does not support C++0x. Minimum required
versions of GCC and Clang are given in the beginning of this document.


Q: "libstdc++.so.6: version `CXXABI_1.3.5' not found" or similar error when running

A: Ensure that the version of libc++ that was used to compile Ragout is similar
to one the you currently using. You can specify an extra search path
to a specific library by setting "LD_LIBRARY_PATH" variable.
