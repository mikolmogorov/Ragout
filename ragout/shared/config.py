#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module stores some configuration parameters
"""

vals =  {
            "overlap" :
            {
                "min_overlap" : 33,
                "max_overlap" : 200,
                "max_path_len" : 30,
                "detect_kmer" : True
            },

            "maf2synteny" :
            [
                (30, 500),
                (100, 5000),
                (500, 50000),
                (5000, 500000)
            ],

            "sibelia" :
            [
                (30, 150),
                (100, 500),
                (500, 1500)
            ],

            "blocks" :
            {
                "small" : [5000, 500, 100],
                "large" : [10000, 500, 100]
            },

            "big_genome_threshold" : 500 * 1024 * 1024,

            "min_synteny_coverage" : 0.6,
            "min_overlap_rate" : 0.5,
            "min_scaffold_gap": 11,
            "max_scaffold_gap": 10000
        }
