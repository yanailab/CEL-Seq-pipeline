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
    """  Run the command, and return the stdout. If there was an error, return zerors. 
         This is because htseq fails on empty files.
    """
    logger.info("ran  : " + " ".join( cmd))
    try:
        out = subprocess.check_output(cmd).splitlines()
    except subprocess.CalledProcessError as e:
        logger.error("HTSeq error with command : %s", e.cmd)
        # for an empty sam file, we want a lot of zeros..
        # but we cannot specify inifinite zeros here because it will take 
        # forever to return them all.
        out = None
    return out


def main(input_files, gff_file, output_dir, extra_params, count_filename, umi="false", procs=50):

    procs = int(procs)
    if umi.lower() in ["true","yes","1"]:
        extra_params += " -u "
    # The first col we need only once, as it is always the same.
    # So we put None, and try to write it on our first chance.
    ht_col1 = None
    ht_col2 = []
    base_names = []
    cmds = []
    for sam_file in input_files:

       base_names += [os.path.splitext(os.path.basename(sam_file))[0]] 
        #  we need base_names for heading the matrix file.
       htseq_cmd =  ["htseq-count-umified"] + shlex.split(extra_params) + [ sam_file, gff_file]
       cmds.append(htseq_cmd)
    
    pool = ThreadPool(procs)
    results = pool.map(run_cmd, cmds)
    for res in results:
        if res is None:
            # HTSeq failed (perhaps empty file), so we put a column of zeros.
            ht_col2.append(cycle(["0"]))
        else:
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
        
