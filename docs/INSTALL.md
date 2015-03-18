Installation instructions for Ragout
====================================

Availability
------------
Ragout was tested under Mac OS and Linux. While it should work
under Windows, we currently do not provide an official support.


Build requirements
------------------
* C++ compiler with C++0x support (GCC 4.6+ / Clang 3.2+ / Apple Clang 4.2+)
* GNU make / Cmake (for building Sibelia)


Runtime depenencies
-------------------

* Python 2.7
* Sibelia [http://github.com/bioinf/Sibelia] or Progresssive Cactus [http://github.com/glennhickey/progressiveCactus]


Binary distribution
-------------------

While we recommend to build Ragout from source on each machine, you also can
use pre-compiled binaries which are available for Linux and Mac OS from 
Releases page on Github: https://github.com/fenderglass/Ragout/releases

In this case, you do not need any installation procedures. If you are unsure,
which version (binary or source) you have downloaded, you can check for
binary files in "lib" directory. If they exist - you have obtained a binary
version, otherwise it is source.


Building
--------

1. To build Ragout native modules, type:
    
        make

2. To build and install Sibelia, use:

        python scripts/install-sibelia.py

If you already have Sibelia installed or do not need it, 
you can skip the second step.


progressiveCactus
-----------------

Ragout can use both Sibelia and Progressive Cactus for synteny block decomposition.
Cactus support is still in early stage, however you already can try it with genomes,
that are too big for Sibelia.

First, download and build Progressive Cactus: https://github.com/glennhickey/progressiveCactus
Then set "CACTUS_INSTALL" environment variable pointing to Cactus installation directory:

    export CACTUS_INSTALL=path_to_cactus


Troubleshooting
---------------

Q: Different errors during compilation, possibly with 
"unrecognized command line option '-std=c++0x'" message:

A: Probably your compiler is too old and does not support C++0x. Minimum
versions of GCC and Clang are mentioned at the beginnig of this documtent.


Q: "libstdc++.so.6: version `CXXABI_1.3.5' not found" or similar error when running

A: Ensure, that version of libc++ which was used to compile Ragout is similar
to one which is used to run it. You can specify an extra search path
to a specific library by setting "LD_LIBRARY_PATH" variable.