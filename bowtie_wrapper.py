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

logger = getLogger('pijp.bowtie_wrapper')

def build_bowtie_command(fastq_file,  index_file, number_of_threads, output_dir, extra_params):


    base_fastq = os.path.splitext(os.path.basename(fastq_file))[0]
    samfile = os.path.join(output_dir, base_fastq + ".sam")

    ##  no-hd means no header lines. 
    ##  -p is for the number of rows.
    bowtie_cmd = "bowtie2  -p {0} {1} {2} -U {3} -S {4} ".format(number_of_threads, extra_params, index_file, fastq_file, samfile)
    return bowtie_cmd


def main(input_files, index_file, number_of_threads, output_dir, bowtie_report_name,extra_params):

    report = []

    filename = os.path.join(output_dir, bowtie_report_name)
    with open(filename, 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(('#sample','total', 'not_aligned', 'aligned_once', 'multi_aligned'))
        for fastq_file in input_files:
            bt2_cmd = build_bowtie_command(fastq_file, index_file, number_of_threads, output_dir, extra_params)
        
            logger.info("ran  :" + bt2_cmd)
            #subprocess.check_call(bt2_cmd, shell=True)
            pro = subprocess.Popen(bt2_cmd, shell=True, stderr = subprocess.PIPE)
            (x, stderr) = pro.communicate()
            assert (pro.returncode == 0 ), "bowtie error %d : %s" % (pro.returncode, stderr)
    
            new_row = ( [fastq_file] + get_stats(stderr.splitlines()) )
         
            writer.writerow(new_row)

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
        total = bt_stderr[1].split()[0]
        not_aligned = bt_stderr[2].split()[0]
        aligned = bt_stderr[3].split()[0]
        multi_aligned =bt_stderr[4].split()[0] 
        return [total, not_aligned, aligned, multi_aligned]
    else:
        return [0,0,0,0]
