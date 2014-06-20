Installation instructions for Ragout
====================================

Availability
------------
Ragout is tested under Mac OS and Linux. While it *should* work
under Windows, we currently do not provide an official support.


Build requirements
------------------
* C++ compiler with C++0x support (GCC 4.6+ / Clang 3.2+ / Apple Clang 4.2+)
* GNU make / Cmake (for building Sibelia)


Runtime depenencies
-------------------

* Python 2.7
* biopython [http://biopython.org]
* Sibelia [https://github.com/bioinf/Sibelia] or Progresssive Cactus [https://github.com/glennhickey/progressiveCactus]


Building
--------

1. To build Ragout native modules type
    
    make

2. To build and install Sibelia, use

    python scripts/install-sibelia.py

If you already have Sibelia installed, you can skip second step.

After this, you can test your installation by typing:

    bin/ragout --help

If you got no errors, installation was successful and you can start using Ragout!


Binary distribution
-------------------

While we recommend to build Ragout from source on each machine, you also can
use pre-compiled binaries which are available for Linux and Mac OS on our
website: <>


progressiveCactus
-----------------

Ragout can use both Sibelia and Progressive Cactus for synteny block decomposition.
Cactus support is still in early stage, however you already can try it with genomes,
that are too bif for Sibelia.

First, download and build Progressive Cactus: https://github.com/glennhickey/progressiveCactus
Then set "CACTUS_INSTALL" environment variable pointing to Cactus installation directory:

    export CACTUS_INSTALL=path_to_cactus


Troubleshooting
---------------

Q: Tons of errors during compilation, possibly with 
"unrecognized command line option '-std=c++0x'" message:

A: Probably your compiler is too old and does not support C++0x. Minimum
versions of GCC and Clang are mentioned at the beginnig of this documtent.