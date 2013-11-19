#!/usr/bin/python2
""" A python wrapper for htseq-count (which is itself in python ofcourse) 
    Should be used via pijpleiding.py

    Outputs one file, the expression matrix.
"""

import subprocess
from logging import getLogger
import shlex
import os.path
import csv
from itertools import cycle

from multiprocessing.pool import ThreadPool

logger = getLogger("pijp.htseq")

def run_cmd(cmd):
    logger.info("ran  : " + " ".join( cmd))
    return subprocess.check_output(cmd).splitlines()

def main(input_files, gff_file, output_dir, extra_params, count_filename, procs=50):

    # The first col we need only once, as it is always the same.
    # So we put None, and try to write it on our first chance.
    ht_col1 = None
    ht_col2 = []
    base_names = []
    cmds = []
    for sam_file in input_files:

       base_names += [os.path.splitext(os.path.basename(sam_file))[0]] 
        #  we need base_names for heading the matrix file.
       htseq_cmd =  ["htseq-count"] + shlex.split(extra_params) + [ sam_file, gff_file]
       cmds.append(htseq_cmd)
    
    pool = ThreadPool(procs)
    results = pool.map(run_cmd, cmds)
    for res in results:
            htout = (x.split('\t') for x in  res)
            htout1, htout2 = zip(*htout)
            if ht_col1 is None:
                ht_col1 = htout1
            ht_col2.append(htout2)
    
    matrix_header = ["#Sample:"] + base_names

    matrix = zip( ht_col1, *ht_col2)
    with open(os.path.join(output_dir, count_filename), "w") as fh:
        matwriter = csv.writer(fh, delimiter='\t')
        matwriter.writerow(matrix_header)
        matwriter.writerows(matrix)
        
