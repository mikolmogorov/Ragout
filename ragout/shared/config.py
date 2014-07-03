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
                "only_on_k" : True
            },

            "maf2synteny" :
            [
                (30, 10),
                (100, 100),
                (500, 1000),
                (1000, 5000),
                (5000, 15000)
            ],

            "sibelia" :
            [
                (30, 150),
                (100, 500),
                (500, 1500)
            ]
        }
