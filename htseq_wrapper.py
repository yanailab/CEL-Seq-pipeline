#!/usr/bin/python2
""" A python wrapper for htseq-count (which is itself in python ofcourse) 
    Should be used via pijpleiding.py

    Outputs one file, the expression matrix.
"""

#import subprocess
from logging import getLogger
import shlex
import os.path
import sys
import csv
from itertools import cycle
import argparse

import htseq_count_umified

from multiprocessing.pool import ThreadPool

def build_argument_opts(cmd_line_params):
    """ Parse the input arguments (string) so that it is easily readable by htseq-count
    """
    #print >> sys.stderr, " ".join(sys.argv)
    parser = argparse.ArgumentParser(description='Python wrapper for htseq-count.')
    # htseq-count arguments
    parser.add_argument('-m','--mode',dest='mode',default='union')
    parser.add_argument('-s','--stranded',dest='stranded',default='yes')
    parser.add_argument('-a','--minaqual',type=int,dest='minaqual',default=0)
    parser.add_argument('-t','--featuretype',dest='featuretype',default='exon')
    parser.add_argument('-i','--idattr',dest='idattr',default='gene_id')
    # more htseq-count commands
    parser.add_argument('-o','--samout',dest='samout',default='')
    parser.add_argument('-q','--quiet',dest='quiet',action="store_true",default=False)
    parser.add_argument('-u', '--umis', dest='umis',action="store_true",default=False)
    #parser.add_argument('-u', '--umis', dest='umis', action="store_true",default='False')
    
    arguments = parser.parse_args(shlex.split(cmd_line_params))
    return arguments

def run_cmd(cmd):
    """  Run the command, and return a feature/count list. 
         If there was an error, return zeros. 
    """
    logger = getLogger("pijp.htseq")
    sam_file = cmd[1]
    gff_file = cmd[2]
    args = build_argument_opts(cmd[0])
    logger.info("ran HTSeq-count: " + sam_file + ', '+ gff_file + ', ' + str(args))
    try:
         out = htseq_count_umified.count_reads_in_features( sam_file, gff_file, args.stranded,
               args.mode, args.featuretype, args.idattr, args.quiet, args.minaqual, 
               args.samout, args.umis)
    except:
         logger.error("HTSeq error with command : %s", cmd)
         out = None
    return out    


def main(input_files, gff_file, output_dir, extra_params, count_filename, umi="false", procs=50):
    """ HTSeq-count wrapper main. Counts input SAMs against given reference with given args.
        Sums all counts in a table.
    """
    procs = int(procs)
    if umi.lower() in ["true","yes","1"]:
        extra_params += " -u "
    # The first col we need only once, as it is always the same.
    # So we put None, and try to write it on our first chance.
    ht_col1 = None
    ht_col2 = []
    base_names = []
    cmds = []     
    # command arguments for each input SAM file
    for sam_file in input_files:
       base_names += [os.path.splitext(os.path.basename(sam_file))[0]] 
        #  we need base_names for heading the matrix file.
       htseq_cmd =  [extra_params, sam_file, gff_file]
       cmds.append(htseq_cmd)
     
    # running htseq-count on multiple processes
    pool = ThreadPool(procs)
    results = pool.map(run_cmd, cmds)
    for res in results:
        if res is None:
            # HTSeq failed (perhaps empty file), so we put a column of zeros.
            ht_col2.append(cycle(["0"]))
        else:
            htout1, htout2 = res
            if ht_col1 is None:
                ht_col1 = htout1
            ht_col2.append(htout2)

    # make a htseq-count matrix of results for output
    matrix_header = ["#Sample:"] + base_names
    matrix = zip( ht_col1, *ht_col2)
    with open(os.path.join(output_dir, count_filename), "w") as fh:
        matwriter = csv.writer(fh, delimiter='\t')
        matwriter.writerow(matrix_header)
        matwriter.writerows(matrix)
        
