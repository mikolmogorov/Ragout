Installation Instructions for Ragout
====================================

Availability
------------
Ragout is available for Mac OS and Linux.


Build Requirements
------------------
* C++ compiler with C++0x support (GCC 4.6+ / Clang 3.2+ / Apple Clang 4.2+)
* GNU make / Cmake (for building Sibelia)


Runtime Depenencies
-------------------

* Python 2.7
* Sibelia [http://github.com/bioinf/Sibelia] or HAL Tools [https://github.com/glennhickey/hal]


Binary Distribution
-------------------

Pre-compiled binaries are available for Linux and Mac OS from 
the releases page [https://github.com/fenderglass/Ragout/releases].
In this case you will not need the installation procedures below.


Building
--------

To build Ragout native modules, type:
    
        make

You will also need either *Sibelia* or *HAL tools* installed (see below)


Sibelia
-------

To build and install Sibelia, use:

        python scripts/install-sibelia.py

If you already have Sieblia installed into your system, it will
be picked up automatically by Ragout.


HAL Tools
---------

HAL alignment produced by *Progressive Cactus* could be used for synteny 
blocks decomposition instead of *Sibelia* (recommended for large genomes). 
If you want to use HAL alignment as input,
you need to install *HAL Tools* package [https://github.com/glennhickey/hal]
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
