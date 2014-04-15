import sys, os, stat
from setuptools import setup, find_packages, Extension
from glob import glob
from platform import uname

compile_args = ["-std=c++11"]
if uname()[0] == "Darwin":                 #a fix for buggy macos
    compile_args.extend(["-stdlib=libc++", "-Qunused-arguments"])

coverlap = Extension("ragout.coverlap",
                    define_macros = [("PYTHON_LIB", 1)],
                    sources = glob("ragout/overlap/cpp_impl/*.cpp"),
                    extra_compile_args = compile_args)

cmaf2synteny = Extension("ragout.cmaf2synteny",
                    define_macros = [("PYTHON_LIB", 1)],
                    sources = glob("ragout/maf2synteny/cpp_impl/*.cpp"),
                    extra_compile_args = compile_args)

#local installation feature
script_args = sys.argv[1:]
INPLACE_EXEC = "ragout_local"
if "inplace" in script_args:
    script_args = ["build_ext", "--inplace"]
    with open(INPLACE_EXEC, "w") as f:
        f.write("#!{0}\nimport ragout.main\nragout.main.main()"
                .format(sys.executable))
    st = os.stat(INPLACE_EXEC)
    os.chmod(INPLACE_EXEC, st.st_mode | stat.S_IEXEC)


setup(
    script_args = script_args,
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
    ext_modules = [coverlap, cmaf2synteny]
)
