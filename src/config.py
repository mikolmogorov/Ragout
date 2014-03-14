import os

RAGOUT_ROOT = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
RAGOUT_LIB_DIR = os.path.join(RAGOUT_ROOT, "lib")

ASSEMBLY_MIN_OVERLAP = 33
ASSEMBLY_MAX_OVERLAP = 100
ASSEMBLY_MAX_PATH_LEN = 6
ASSEMBLY_KMER = 54

CACTUS_DIR = "/home/volrath/Bioinf/tools/progressiveCactus/"
