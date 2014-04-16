import sys, os, stat
from setuptools import setup, find_packages, Extension
from glob import glob
from platform import uname


#extensions
compile_args = ["-std=c++0x"]
if uname()[0] == "Darwin":                 #a fix for buggy macos
    compile_args.extend(["-stdlib=libc++", "-Qunused-arguments"])


coverlap = Extension("ragout.coverlap",
                    language = "c++",
                    define_macros = [("PYTHON_LIB", 1)],
                    sources = glob("ragout/overlap/cpp_impl/*.cpp"),
                    extra_compile_args = compile_args)

cmaf2synteny = Extension("ragout.cmaf2synteny",
                    language = "c++",
                    define_macros = [("PYTHON_LIB", 1)],
                    sources = glob("ragout/maf2synteny/cpp_impl/*.cpp"),
                    extra_compile_args = compile_args)


#local installation feature
script_args = sys.argv[1:]
INPLACE_EXEC = "run-ragout"
if "inplace" in script_args:
    script_args = ["build_ext", "--inplace"]
    open(INPLACE_EXEC, "w").write("#!" + sys.executable + "\n"
                                  "from ragout.main import main\n"
                                  "main()")
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
            "run-ragout = ragout.main:main"
        ]
    },
    install_requires = ["biopython", "networkx>=1.8"],
    ext_modules = [coverlap, cmaf2synteny],
    data_files = [("share/ragout/doc", glob("doc/*"))]
                  #("share/ragout/examples", get_files("examples"))]
)
