from setuptools import setup, find_packages, Extension
from glob import glob
from platform import uname

compile_args = ["-std=c++11"]
if uname()[0] == "Darwin":                 #a fix for buggy macos
    compile_args.extend(["-stdlib=libc++", "-Qunused-arguments"])

coverlap = Extension("coverlap",
                    define_macros = [("PYTHON_LIB", 1)],
                    sources = glob("ragout/overlap/cpp_impl/*.cpp"),
                    extra_compile_args = compile_args)

cmaf2synteny = Extension("cmaf2synteny",
                    define_macros = [("PYTHON_LIB", 1)],
                    sources = glob("ragout/maf2synteny/cpp_impl/*.cpp"),
                    extra_compile_args = compile_args)

setup(
    name = "ragout",
    version = "0.2b",
    author = "Mikhail Kolmogorov",
    author_email = "fenderglass@gmail.com",
    description = "A tool for reference-assisted assembly",
    license = "GPLv2",
    keywords = "bioinformatics",
    url = "http://github.com/fenderglass/Ragout",
    packages = find_packages(),
    long_description = open("README.md", "r").read(),
    entry_points = {
        "console_scripts" : [
            "ragout = ragout.main:main",
            "maf2synteny = ragout.maf2synteny.maf2synteny:main"
        ]
    },
    install_requires = ["biopython", "networkx"],
    ext_modules = [coverlap, cmaf2synteny],
    #package_data = {"" : glob("doc/*")},
    #scripts = glob("scripts/*.py"),
    #include_package_data=True
    #classifiers=[
    #    "Development Status :: 3 - Alpha",
    #    "Topic :: Utilities",
    #    "License :: OSI Approved :: BSD License",
    #],
)
