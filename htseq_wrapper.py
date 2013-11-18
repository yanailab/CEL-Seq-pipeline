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

logger = getLogger("pijp.htseq")

def main(input_files, gff_file, output_dir, extra_params, count_filename):

    # The first col we need only once, as it is always the same.
    # So we put None, and try to write it on our first chance.
    ht_col1 = None
    ht_col2 = []
    base_names = []
    processes = []
    for sam_file in input_files:

        #  we need base_names for heading the matrix file.
        base_names += [os.path.splitext(os.path.basename(sam_file))[0]] 

        htseq_cmd =  ["htseq-count"] + shlex.split(extra_params) + [ sam_file, gff_file]
        logger.info("ran  : " + " ".join( htseq_cmd))
        processes.append( subprocess.Popen((htseq_cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE))

    for pro in processes:
        (stdout, stderr) = pro.communicate()
        if pro.returncode == 0:
            htout = (x.split('\t') for x in  stdout.splitlines())
            htout1, htout2 = zip(*htout)
            if ht_col1 is None:
                ht_col1 = htout1
            ht_col2.append(htout2)
        else:
            if os.stat(sam_file).st_size == 0:
                logger.error(" Empty sam file : %s ", sam_file)
                # for an empty sam file, we want a lot of zeros..
                ht_col2.append(cycle("0"))
            else:
                raise  # some other error we don't know about.
    
    matrix_header = ["#Sample:"] + base_names

    matrix = zip( ht_col1, *ht_col2)
    with open(os.path.join(output_dir, count_filename), "w") as fh:
        matwriter = csv.writer(fh, delimiter='\t')
        matwriter.writerow(matrix_header)
        matwriter.writerows(matrix)
        
