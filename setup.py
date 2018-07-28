#(c) 2016 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

from __future__ import print_function
import sys

#Check Python version
if sys.version_info[:2] != (2, 7):
    print("Error: Flye requires Python version 2.7 ({0}.{1} detected)."
          .format(sys.version_info[0], sys.version_info[1]))
    sys.exit(-1)

from distutils.core import setup
from distutils.command.build import build as DistutilsBuild
import subprocess

from ragout.__version__ import __version__


class MakeBuild(DistutilsBuild):
    def run(self):
        try:
            subprocess.check_call(['make'])
        except subprocess.CalledProcessError as e:
            print ("Compilation error: ", e)
            return
        DistutilsBuild.run(self)

setup(name='ragout',
      version=__version__,
      description='A tool for chromosome assembly using multiple references',
      url='https://github.com/fenderglass/Ragout',
      author='Mikhail Kolmogorov',
      author_email = 'fenderglass@gmail.com',
      license='BSD-3-Clause',
      packages=['ragout', 'ragout/assembly_graph', 'ragout/breakpoint_graph',
                'ragout/maf2synteny', 'ragout/overlap', 'ragout/parsers',
                'ragout/phylogeny', 'ragout/scaffolder', 'ragout/shared',
                'ragout/synteny_backend'],
      package_data={'ragout': ['LICENSE']},
      scripts = ['bin/ragout-maf2synteny', 'bin/ragout-overlap', 'bin/ragout'],
      cmdclass={'build': MakeBuild}
      )
