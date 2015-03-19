Installation Instructions for Ragout
====================================

Availability
------------
Ragout was tested under Mac OS and Linux. While it should work
under Windows, we currently do not provide an official support.


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


HAL Tools
---------

HAL alignment could be used for synteny blocks decomposition instead of Sibelia
(recommended for large genomes). If you want to use HAL alignment as input,
you need to install HAL Tools package: https://github.com/glennhickey/hal.
Follow the manuals and do not forget to properly set PATH and PYTHONPATH
environment variables as it is described in them.


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