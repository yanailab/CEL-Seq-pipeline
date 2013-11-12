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

logger = getLogger("pijp.htseq")

def main(input_files, gff_file, output_dir, extra_params, count_filename):

    htout = []
    base_names = []
    for sam_file in input_files:

        #  we need base_names for heading the matrix file.
        base_names += [os.path.splitext(os.path.basename(sam_file))[0]] 

        htseq_cmd =  ["htseq-count"] + shlex.split(extra_params) + [ sam_file, gff_file]
        logger.info("ran  : " + " ".join( htseq_cmd))
        try:
            htout += [subprocess.check_output(htseq_cmd).splitlines()]
        except subprocess.CalledProcessError:
            if os.stat(sam_file).st_size == 0:
                logger.error(" Empty sam file : %s ", sam_file)
            else:
                raise
    
    matrix = [["#Sample:"]+base_names]
    for pairs in zip(*htout):
        first_col = pairs[0].split('\t')[0]  # first half of first pair
        second_halves = [pair.split('\t')[1] for pair in pairs]
        matrix += [[first_col] + second_halves] 
    with open(os.path.join(output_dir, count_filename), "w") as fh:
        matwriter = csv.writer(fh, delimiter='\t')
        matwriter.writerows(matrix)
        
