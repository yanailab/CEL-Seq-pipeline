#!/usr/bin/python2
""" run bowtie with specified parameter file 
"""

import subprocess
from logging import getLogger
import os
import shlex
import argparse
import csv
from glob import glob
from multiprocessing.pool import ThreadPool

logger = getLogger('pijp.bowtie_wrapper')

def build_bowtie_command(fastq_file,  index_file, number_of_threads, output_dir, extra_params):
    base_fastq = os.path.splitext(os.path.basename(fastq_file))[0]
    samfile = os.path.join(output_dir, base_fastq + ".sam")

    ##  no-hd means no header lines. 
    ##  -p is for the number of rows.
    bowtie_cmd = "bowtie2  -p {0} {1} {2} -U {3} -S {4} ".format(number_of_threads, extra_params, index_file, fastq_file, samfile)
    return bowtie_cmd

def run_cmd((cmd,fastq_file)):
    """  Run the command, and return the stderr. 
    """
    logger.info("ran  : " + cmd)
    pro = subprocess.Popen(cmd, shell=True, stderr = subprocess.PIPE)
    (x, stderr) = pro.communicate()
    assert (pro.returncode == 0 ), "bowtie error %d : %s" % (pro.returncode, stderr)
    new_row = ( [fastq_file] + get_stats(stderr.splitlines()) )
    logger.info("finished  : " + cmd)
    return new_row

def main(input_files, index_file, number_of_threads, output_dir, bowtie_report_name,extra_params,procs=10):
    report = []
    ht_col2 = []
    base_names = []
    cmds = []
    for fastq_file in input_files:
       base_names += [os.path.splitext(os.path.basename(fastq_file))[0]]
       #  we need base_names for heading the matrix file.
       bt2_cmd = build_bowtie_command(fastq_file, index_file, number_of_threads, output_dir, extra_params)       
       cmds.append((bt2_cmd,fastq_file))
    ht_col1 = "\n".join(base_names)
    pool = ThreadPool(int(procs))
    results = pool.map(run_cmd, cmds)
    for res in results:
        htout = "\t".join(res)
        ht_col2.append(htout)
    matrix_header = "\t".join(['#sample','total', 'not_aligned', 'aligned_once', 'multi_aligned', '% mapped']) + "\n"
    matrix = "\n".join(ht_col2)
    filename = os.path.join(output_dir, bowtie_report_name)
    with open(filename, 'wb') as f:
        f.write(matrix_header)
        f.write(matrix)

def get_stats(bt_stderr):
    ## Write all warnings to log, and skip them.
    # For millions of warnings this takes too long..
    #while bt_stderr[0].startswith("Warning"):
    #    logger.info("BOWTIE : " + bt_stderr[0])
    #    bt_stderr.pop(0)
    if len(bt_stderr)  > 5:
        logger.info("Bowtie probably had {0} warnings [e.g. read too small]".format(len(bt_stderr)-4))
        bt_stderr = bt_stderr[-5:]
    if int(bt_stderr[0].split()[0]) != 0 :
        total = str(int(bt_stderr[0].split()[0]))
        not_aligned = str(int(bt_stderr[1].split()[0]))
        aligned = str(int(bt_stderr[2].split()[0]))
        multi_aligned = str(int(bt_stderr[3].split()[0]))
        mapped_percent = bt_stderr[4].split()[0]
        return [total, not_aligned, aligned, multi_aligned, mapped_percent]
    else:
        return ['0','0','0','0','0']
